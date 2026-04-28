//! History tracking for topology operations, matching OCCT BRepAlgoAPI_BuilderShape capabilities.
//!
//! This module provides comprehensive tracking of:
//! - **Modifications**: Entities modified during an operation
//! - **Generations**: New entities created during an operation
//! - **Deletions**: Entities removed during an operation
//! - **Queries**: Ancestor/descendant relationships across operations
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::history::{BooleanHistory, HistoryTracker, GenerationCause, DeletionReason};
//!
//! // Create a history tracker for a boolean operation
//! let mut tracker = HistoryTracker::new();
//!
//! // Record modifications, generations, and deletions
//! tracker.record_face_modified(0, 1);  // Face 0 modified to face 1
//! tracker.record_face_generated(2, GenerationCause::Intersection);     // Face 2 was newly generated
//! tracker.record_face_deleted(3, DeletionReason::Custom("Removed by boolean cut".to_string()));
//!
//! // Query the history
//! assert!(tracker.has_modified());
//! assert!(tracker.has_generated());
//! assert!(tracker.has_deleted());
//! ```

use rcad_kernel::{
    persistent_naming::{
        NamingEvent, OperationStats, PersistentId, PersistentNamingEngine,
    },
    BRep, PersistentNamingHooks, TopoEntityRef,
};
use std::collections::HashMap;

/// Types of boolean operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOperationType {
    Union,
    Intersection,
    Difference,
}

impl From<BooleanOperationType> for rcad_kernel::persistent_naming::OperationType {
    fn from(op: BooleanOperationType) -> Self {
        match op {
            BooleanOperationType::Union => rcad_kernel::persistent_naming::OperationType::BooleanUnion,
            BooleanOperationType::Intersection => rcad_kernel::persistent_naming::OperationType::BooleanIntersection,
            BooleanOperationType::Difference => rcad_kernel::persistent_naming::OperationType::BooleanDifference,
        }
    }
}

/// Reason for entity deletion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DeletionReason {
    /// Removed by boolean operation.
    BooleanOperation,
    /// Removed because it was outside the result volume.
    OutsideResult,
    /// Removed due to overlapping geometry.
    Overlap,
    /// Removed due to tolerance issues.
    Tolerance,
    /// Removed during healing/repair.
    Healing,
    /// Custom reason with description.
    Custom(String),
}

impl DeletionReason {
    /// Returns a human-readable description of the deletion reason.
    pub fn description(&self) -> &str {
        match self {
            DeletionReason::BooleanOperation => "Removed by boolean operation",
            DeletionReason::OutsideResult => "Outside result volume",
            DeletionReason::Overlap => "Overlapping geometry",
            DeletionReason::Tolerance => "Tolerance issues",
            DeletionReason::Healing => "Removed during healing",
            DeletionReason::Custom(s) => s.as_str(),
        }
    }
}

/// Record of a deleted entity.
#[derive(Debug, Clone)]
pub struct DeletionRecord {
    /// Index of the deleted entity.
    pub entity_index: usize,
    /// Type of the deleted entity.
    pub entity_type: EntityType,
    /// Reason for deletion.
    pub reason: DeletionReason,
    /// Source input (A or B) if applicable.
    pub source: Option<InputSource>,
}

/// Types of topological entities.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EntityType {
    Vertex,
    Edge,
    Face,
    Shell,
    Solid,
}

/// Source of an entity in a boolean operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InputSource {
    /// Entity came from input A.
    A,
    /// Entity came from input B.
    B,
    /// Entity was generated during the operation.
    Generated,
}

/// Tracks the origin of each face in a boolean operation result.
///
/// Analogous to OCCT `BRepAlgoAPI_BuilderShape::Modified/Generated/Deleted`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceOrigin {
    /// Face came from solid A; value is the DS face index of the source face.
    FromA(usize),
    /// Face came from solid B; value is the DS face index of the source face.
    FromB(usize),
    /// Face was generated at the intersection boundary (not yet produced by this
    /// implementation — reserved for future use).
    Generated,
}

/// Tracks the origin of each edge in a boolean operation result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeOrigin {
    /// Edge came from solid A; value is the original edge index in solid A.
    FromA(usize),
    /// Edge came from solid B; value is the original edge index in solid B.
    FromB(usize),
    /// Edge is an intersection edge generated at the boolean boundary.
    Generated,
    /// Edge was created from a partial (split) segment of an original edge in A.
    SplitFromA(usize),
    /// Edge was created from a partial (split) segment of an original edge in B.
    SplitFromB(usize),
}

/// Tracks the origin of each vertex in a boolean operation result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VertexOrigin {
    /// Vertex came from solid A; value is the original vertex index in solid A.
    FromA(usize),
    /// Vertex came from solid B; value is the original vertex index in solid B.
    FromB(usize),
    /// Vertex was created at an A-B intersection point.
    Intersection,
}

/// Aggregate origin of a result shell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellOrigin {
    /// Every tracked face in the shell came from solid A.
    FromA,
    /// Every tracked face in the shell came from solid B.
    FromB,
    /// Every tracked face in the shell was generated.
    Generated,
    /// The shell contains a mixture of A/B/generated faces.
    Mixed,
}

/// Aggregate origin of a result solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidOrigin {
    /// Every tracked shell in the solid came from solid A.
    FromA,
    /// Every tracked shell in the solid came from solid B.
    FromB,
    /// Every tracked shell in the solid was generated.
    Generated,
    /// The solid contains a mixture of A/B/generated shells or faces.
    Mixed,
}

/// Record of an entity modification.
#[derive(Debug, Clone)]
pub struct ModificationRecord {
    /// Source entity index.
    pub source_index: usize,
    /// Source input (A or B).
    pub source: InputSource,
    /// Result entity indices (may be multiple for splits).
    pub result_indices: Vec<usize>,
    /// Type of modification.
    pub modification_type: ModificationType,
}

/// Types of modifications that can occur.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModificationType {
    /// Entity was directly preserved (same geometry, possibly new index).
    Preserved,
    /// Entity was split into multiple entities.
    Split,
    /// Entity was merged with others.
    Merged,
    /// Entity geometry was modified.
    GeometryModified,
}

/// Record of a generated entity.
#[derive(Debug, Clone)]
pub struct GenerationRecord {
    /// Index of the generated entity.
    pub entity_index: usize,
    /// Type of the generated entity.
    pub entity_type: EntityType,
    /// What caused this entity to be generated.
    pub cause: GenerationCause,
    /// Parent entity indices that contributed to generation (if any).
    pub parent_indices: Vec<usize>,
}

/// Causes for entity generation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GenerationCause {
    /// Generated at an intersection.
    Intersection,
    /// Generated to fill a gap.
    GapFill,
    /// Generated as a new boundary.
    NewBoundary,
    /// Generated during healing.
    HealingRepair,
    /// Custom reason.
    Custom(String),
}

// ─────────────────────────────────────────────────────────────────────────────
// HistoryTracker
// ─────────────────────────────────────────────────────────────────────────────

/// Comprehensive tracker for topology operation history.
///
/// This struct provides OCCT BRepAlgoAPI_BuilderShape-like capabilities for:
/// - Tracking modifications (`Modified`, `IsModified`, `HasModified`)
/// - Tracking generations (`Generated`, `HasGenerated`)
/// - Tracking deletions (`IsDeleted`, `HasDeleted`)
/// - Ancestor/descendant queries across operations
#[derive(Debug, Clone, Default)]
pub struct HistoryTracker {
    /// Face modification records: source_index -> result_indices.
    face_modifications: HashMap<usize, ModificationRecord>,
    /// Edge modification records: source_index -> result_indices.
    edge_modifications: HashMap<usize, ModificationRecord>,
    /// Vertex modification records: source_index -> result_indices.
    vertex_modifications: HashMap<usize, ModificationRecord>,
    /// Generated faces with their generation cause.
    generated_faces: Vec<GenerationRecord>,
    /// Generated edges with their generation cause.
    generated_edges: Vec<GenerationRecord>,
    /// Generated vertices with their generation cause.
    generated_vertices: Vec<GenerationRecord>,
    /// Deleted entities with reasons.
    deleted_entities: Vec<DeletionRecord>,
    /// Quick lookup for deleted status: (entity_type, index) -> deletion_index.
    deleted_lookup: HashMap<(EntityType, usize), usize>,
    /// Reverse lookup: result entity -> source entity.
    result_to_source: HashMap<(EntityType, usize), (InputSource, usize)>,
    /// Source input tracking for each entity type.
    face_sources: HashMap<usize, (InputSource, usize)>,
    edge_sources: HashMap<usize, (InputSource, usize)>,
    vertex_sources: HashMap<usize, (InputSource, usize)>,
}

impl HistoryTracker {
    /// Create a new empty history tracker.
    pub fn new() -> Self {
        Self::default()
    }

    // ── Modification Recording ───────────────────────────────────────────────────

    /// Record that a face was modified from a source to one or more results.
    pub fn record_face_modified(&mut self, source_index: usize, result_index: usize) {
        self.record_face_modified_multi(source_index, vec![result_index], InputSource::A);
    }

    /// Record that a face was modified from a source to multiple results (split).
    pub fn record_face_modified_multi(
        &mut self,
        source_index: usize,
        result_indices: Vec<usize>,
        source: InputSource,
    ) {
        let mod_type = if result_indices.len() > 1 {
            ModificationType::Split
        } else {
            ModificationType::Preserved
        };

        let record = ModificationRecord {
            source_index,
            source,
            result_indices: result_indices.clone(),
            modification_type: mod_type,
        };

        self.face_modifications.insert(source_index, record);

        // Update reverse lookup.
        for &result_idx in &result_indices {
            self.face_sources.insert(result_idx, (source, source_index));
            self.result_to_source.insert((EntityType::Face, result_idx), (source, source_index));
        }
    }

    /// Record that an edge was modified from a source to one or more results.
    pub fn record_edge_modified(
        &mut self,
        source_index: usize,
        result_indices: Vec<usize>,
        source: InputSource,
        mod_type: ModificationType,
    ) {
        let record = ModificationRecord {
            source_index,
            source,
            result_indices: result_indices.clone(),
            modification_type: mod_type,
        };

        self.edge_modifications.insert(source_index, record);

        // Update reverse lookup.
        for &result_idx in &result_indices {
            self.edge_sources.insert(result_idx, (source, source_index));
            self.result_to_source.insert((EntityType::Edge, result_idx), (source, source_index));
        }
    }

    /// Record that a vertex was modified from a source to one or more results.
    pub fn record_vertex_modified(
        &mut self,
        source_index: usize,
        result_indices: Vec<usize>,
        source: InputSource,
    ) {
        let mod_type = if result_indices.len() > 1 {
            ModificationType::Split
        } else {
            ModificationType::Preserved
        };

        let record = ModificationRecord {
            source_index,
            source,
            result_indices: result_indices.clone(),
            modification_type: mod_type,
        };

        self.vertex_modifications.insert(source_index, record);

        // Update reverse lookup.
        for &result_idx in &result_indices {
            self.vertex_sources.insert(result_idx, (source, source_index));
            self.result_to_source.insert((EntityType::Vertex, result_idx), (source, source_index));
        }
    }

    // ── Generation Recording ─────────────────────────────────────────────────────

    /// Record a generated face.
    pub fn record_face_generated(&mut self, face_index: usize, cause: GenerationCause) {
        self.record_face_generated_with_parents(face_index, cause, vec![]);
    }

    /// Record a generated face with parent indices.
    pub fn record_face_generated_with_parents(
        &mut self,
        face_index: usize,
        cause: GenerationCause,
        parent_indices: Vec<usize>,
    ) {
        let record = GenerationRecord {
            entity_index: face_index,
            entity_type: EntityType::Face,
            cause,
            parent_indices,
        };
        self.generated_faces.push(record);
        self.face_sources.insert(face_index, (InputSource::Generated, face_index));
    }

    /// Record a generated edge (e.g., intersection edge).
    pub fn record_edge_generated(&mut self, edge_index: usize, cause: GenerationCause) {
        self.record_edge_generated_with_parents(edge_index, cause, vec![]);
    }

    /// Record a generated edge with parent indices.
    pub fn record_edge_generated_with_parents(
        &mut self,
        edge_index: usize,
        cause: GenerationCause,
        parent_indices: Vec<usize>,
    ) {
        let record = GenerationRecord {
            entity_index: edge_index,
            entity_type: EntityType::Edge,
            cause,
            parent_indices,
        };
        self.generated_edges.push(record);
        self.edge_sources.insert(edge_index, (InputSource::Generated, edge_index));
    }

    /// Record a generated vertex (e.g., intersection vertex).
    pub fn record_vertex_generated(&mut self, vertex_index: usize, cause: GenerationCause) {
        self.record_vertex_generated_with_parents(vertex_index, cause, vec![]);
    }

    /// Record a generated vertex with parent indices.
    pub fn record_vertex_generated_with_parents(
        &mut self,
        vertex_index: usize,
        cause: GenerationCause,
        parent_indices: Vec<usize>,
    ) {
        let record = GenerationRecord {
            entity_index: vertex_index,
            entity_type: EntityType::Vertex,
            cause,
            parent_indices,
        };
        self.generated_vertices.push(record);
        self.vertex_sources.insert(vertex_index, (InputSource::Generated, vertex_index));
    }

    // ── Deletion Recording ───────────────────────────────────────────────────────

    /// Record a deleted face.
    pub fn record_face_deleted(&mut self, face_index: usize, reason: DeletionReason) {
        self.record_entity_deleted(face_index, EntityType::Face, reason, None);
    }

    /// Record a deleted face with source info.
    pub fn record_face_deleted_with_source(
        &mut self,
        face_index: usize,
        reason: DeletionReason,
        source: InputSource,
    ) {
        self.record_entity_deleted(face_index, EntityType::Face, reason, Some(source));
    }

    /// Record a deleted edge.
    pub fn record_edge_deleted(&mut self, edge_index: usize, reason: DeletionReason) {
        self.record_entity_deleted(edge_index, EntityType::Edge, reason, None);
    }

    /// Record a deleted edge with source info.
    pub fn record_edge_deleted_with_source(
        &mut self,
        edge_index: usize,
        reason: DeletionReason,
        source: InputSource,
    ) {
        self.record_entity_deleted(edge_index, EntityType::Edge, reason, Some(source));
    }

    /// Record a deleted vertex.
    pub fn record_vertex_deleted(&mut self, vertex_index: usize, reason: DeletionReason) {
        self.record_entity_deleted(vertex_index, EntityType::Vertex, reason, None);
    }

    /// Record a deleted vertex with source info.
    pub fn record_vertex_deleted_with_source(
        &mut self,
        vertex_index: usize,
        reason: DeletionReason,
        source: InputSource,
    ) {
        self.record_entity_deleted(vertex_index, EntityType::Vertex, reason, Some(source));
    }

    /// Record a deleted entity.
    pub fn record_entity_deleted(
        &mut self,
        entity_index: usize,
        entity_type: EntityType,
        reason: DeletionReason,
        source: Option<InputSource>,
    ) {
        let deletion_index = self.deleted_entities.len();
        let record = DeletionRecord {
            entity_index,
            entity_type,
            reason,
            source,
        };
        self.deleted_entities.push(record);
        self.deleted_lookup.insert((entity_type, entity_index), deletion_index);
    }

    // ── Modification Queries (OCCT-style) ────────────────────────────────────────

    /// Returns true if any entities were modified during the operation.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::HasModified()`.
    pub fn has_modified(&self) -> bool {
        !self.face_modifications.is_empty()
            || !self.edge_modifications.is_empty()
            || !self.vertex_modifications.is_empty()
    }

    /// Returns true if a specific face was modified.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::IsModified()`.
    pub fn is_face_modified(&self, source_index: usize) -> bool {
        self.face_modifications.contains_key(&source_index)
    }

    /// Returns true if a specific edge was modified.
    pub fn is_edge_modified(&self, source_index: usize) -> bool {
        self.edge_modifications.contains_key(&source_index)
    }

    /// Returns true if a specific vertex was modified.
    pub fn is_vertex_modified(&self, source_index: usize) -> bool {
        self.vertex_modifications.contains_key(&source_index)
    }

    /// Returns the list of result faces that came from a source face.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::Modified()`.
    pub fn modified_faces(&self, source_index: usize) -> Vec<usize> {
        self.face_modifications
            .get(&source_index)
            .map(|r| r.result_indices.clone())
            .unwrap_or_default()
    }

    /// Returns the list of result edges that came from a source edge.
    pub fn modified_edges(&self, source_index: usize) -> Vec<usize> {
        self.edge_modifications
            .get(&source_index)
            .map(|r| r.result_indices.clone())
            .unwrap_or_default()
    }

    /// Returns the list of result vertices that came from a source vertex.
    pub fn modified_vertices(&self, source_index: usize) -> Vec<usize> {
        self.vertex_modifications
            .get(&source_index)
            .map(|r| r.result_indices.clone())
            .unwrap_or_default()
    }

    /// Returns the modification record for a face, if any.
    pub fn face_modification_record(&self, source_index: usize) -> Option<&ModificationRecord> {
        self.face_modifications.get(&source_index)
    }

    /// Returns the modification record for an edge, if any.
    pub fn edge_modification_record(&self, source_index: usize) -> Option<&ModificationRecord> {
        self.edge_modifications.get(&source_index)
    }

    /// Returns the modification record for a vertex, if any.
    pub fn vertex_modification_record(&self, source_index: usize) -> Option<&ModificationRecord> {
        self.vertex_modifications.get(&source_index)
    }

    // ── Generation Queries (OCCT-style) ───────────────────────────────────────────

    /// Returns true if any entities were generated during the operation.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::HasGenerated()`.
    pub fn has_generated(&self) -> bool {
        !self.generated_faces.is_empty()
            || !self.generated_edges.is_empty()
            || !self.generated_vertices.is_empty()
    }

    /// Returns true if a specific result face was generated.
    pub fn is_face_generated(&self, result_index: usize) -> bool {
        self.generated_faces.iter().any(|r| r.entity_index == result_index)
    }

    /// Returns true if a specific result edge was generated.
    pub fn is_edge_generated(&self, result_index: usize) -> bool {
        self.generated_edges.iter().any(|r| r.entity_index == result_index)
    }

    /// Returns true if a specific result vertex was generated.
    pub fn is_vertex_generated(&self, result_index: usize) -> bool {
        self.generated_vertices.iter().any(|r| r.entity_index == result_index)
    }

    /// Returns all generated faces.
    pub fn generated_faces(&self) -> &[GenerationRecord] {
        &self.generated_faces
    }

    /// Returns all generated edges.
    pub fn generated_edges(&self) -> &[GenerationRecord] {
        &self.generated_edges
    }

    /// Returns all generated vertices.
    pub fn generated_vertices(&self) -> &[GenerationRecord] {
        &self.generated_vertices
    }

    /// Returns the generation record for a face, if it was generated.
    pub fn face_generation_record(&self, face_index: usize) -> Option<&GenerationRecord> {
        self.generated_faces.iter().find(|r| r.entity_index == face_index)
    }

    /// Returns the generation record for an edge, if it was generated.
    pub fn edge_generation_record(&self, edge_index: usize) -> Option<&GenerationRecord> {
        self.generated_edges.iter().find(|r| r.entity_index == edge_index)
    }

    /// Returns the generation record for a vertex, if it was generated.
    pub fn vertex_generation_record(&self, vertex_index: usize) -> Option<&GenerationRecord> {
        self.generated_vertices.iter().find(|r| r.entity_index == vertex_index)
    }

    // ── Deletion Queries (OCCT-style) ─────────────────────────────────────────────

    /// Returns true if any entities were deleted during the operation.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::HasDeleted()`.
    pub fn has_deleted(&self) -> bool {
        !self.deleted_entities.is_empty()
    }

    /// Returns true if a specific entity was deleted.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::IsDeleted()`.
    pub fn is_deleted(&self, entity_index: usize, entity_type: EntityType) -> bool {
        self.deleted_lookup.contains_key(&(entity_type, entity_index))
    }

    /// Returns true if a specific face was deleted.
    pub fn is_face_deleted(&self, face_index: usize) -> bool {
        self.is_deleted(face_index, EntityType::Face)
    }

    /// Returns true if a specific edge was deleted.
    pub fn is_edge_deleted(&self, edge_index: usize) -> bool {
        self.is_deleted(edge_index, EntityType::Edge)
    }

    /// Returns true if a specific vertex was deleted.
    pub fn is_vertex_deleted(&self, vertex_index: usize) -> bool {
        self.is_deleted(vertex_index, EntityType::Vertex)
    }

    /// Returns the deletion record for an entity, if it was deleted.
    pub fn deletion_record(&self, entity_index: usize, entity_type: EntityType) -> Option<&DeletionRecord> {
        self.deleted_lookup
            .get(&(entity_type, entity_index))
            .map(|&idx| &self.deleted_entities[idx])
    }

    /// Returns all deleted entities.
    pub fn deleted_entities(&self) -> &[DeletionRecord] {
        &self.deleted_entities
    }

    /// Returns deleted faces.
    pub fn deleted_faces(&self) -> impl Iterator<Item = &DeletionRecord> {
        self.deleted_entities.iter().filter(|r| r.entity_type == EntityType::Face)
    }

    /// Returns deleted edges.
    pub fn deleted_edges(&self) -> impl Iterator<Item = &DeletionRecord> {
        self.deleted_entities.iter().filter(|r| r.entity_type == EntityType::Edge)
    }

    /// Returns deleted vertices.
    pub fn deleted_vertices(&self) -> impl Iterator<Item = &DeletionRecord> {
        self.deleted_entities.iter().filter(|r| r.entity_type == EntityType::Vertex)
    }

    // ── Ancestor/Descendant Queries ──────────────────────────────────────────────

    /// Get the source (original) entity for a result entity.
    /// Returns (source, source_index) if found.
    pub fn get_source(&self, entity_type: EntityType, result_index: usize) -> Option<(InputSource, usize)> {
        match entity_type {
            EntityType::Face => self.face_sources.get(&result_index).copied(),
            EntityType::Edge => self.edge_sources.get(&result_index).copied(),
            EntityType::Vertex => self.vertex_sources.get(&result_index).copied(),
            _ => None,
        }
    }

    /// Get all result entities derived from a source entity.
    pub fn get_results(&self, entity_type: EntityType, source_index: usize) -> Vec<usize> {
        match entity_type {
            EntityType::Face => self.modified_faces(source_index),
            EntityType::Edge => self.modified_edges(source_index),
            EntityType::Vertex => self.modified_vertices(source_index),
            _ => vec![],
        }
    }

    /// Get the input source (A, B, or Generated) for an entity.
    pub fn get_input_source(&self, entity_type: EntityType, index: usize) -> Option<InputSource> {
        match entity_type {
            EntityType::Face => self.face_sources.get(&index).map(|(s, _)| *s),
            EntityType::Edge => self.edge_sources.get(&index).map(|(s, _)| *s),
            EntityType::Vertex => self.vertex_sources.get(&index).map(|(s, _)| *s),
            _ => None,
        }
    }

    /// Count entities by input source.
    pub fn count_by_source(&self, entity_type: EntityType, source: InputSource) -> usize {
        let sources = match entity_type {
            EntityType::Face => &self.face_sources,
            EntityType::Edge => &self.edge_sources,
            EntityType::Vertex => &self.vertex_sources,
            _ => return 0,
        };
        sources.values().filter(|(s, _)| *s == source).count()
    }

    /// Get statistics about the tracked history.
    pub fn statistics(&self) -> HistoryStatistics {
        HistoryStatistics {
            modified_faces: self.face_modifications.len(),
            modified_edges: self.edge_modifications.len(),
            modified_vertices: self.vertex_modifications.len(),
            generated_faces: self.generated_faces.len(),
            generated_edges: self.generated_edges.len(),
            generated_vertices: self.generated_vertices.len(),
            deleted_faces: self.deleted_faces().count(),
            deleted_edges: self.deleted_edges().count(),
            deleted_vertices: self.deleted_vertices().count(),
        }
    }

    /// Clear all tracking data.
    pub fn clear(&mut self) {
        self.face_modifications.clear();
        self.edge_modifications.clear();
        self.vertex_modifications.clear();
        self.generated_faces.clear();
        self.generated_edges.clear();
        self.generated_vertices.clear();
        self.deleted_entities.clear();
        self.deleted_lookup.clear();
        self.result_to_source.clear();
        self.face_sources.clear();
        self.edge_sources.clear();
        self.vertex_sources.clear();
    }

    /// Merge another history tracker into this one.
    /// Useful for combining histories from multiple operations.
    pub fn merge(&mut self, other: &HistoryTracker) {
        // Merge modifications.
        for (idx, record) in &other.face_modifications {
            self.face_modifications.insert(*idx, record.clone());
        }
        for (idx, record) in &other.edge_modifications {
            self.edge_modifications.insert(*idx, record.clone());
        }
        for (idx, record) in &other.vertex_modifications {
            self.vertex_modifications.insert(*idx, record.clone());
        }

        // Merge generations.
        self.generated_faces.extend(other.generated_faces.clone());
        self.generated_edges.extend(other.generated_edges.clone());
        self.generated_vertices.extend(other.generated_vertices.clone());

        // Merge deletions.
        for record in &other.deleted_entities {
            self.record_entity_deleted(
                record.entity_index,
                record.entity_type,
                record.reason.clone(),
                record.source,
            );
        }

        // Merge source lookups.
        for (idx, src) in &other.face_sources {
            self.face_sources.insert(*idx, *src);
        }
        for (idx, src) in &other.edge_sources {
            self.edge_sources.insert(*idx, *src);
        }
        for (idx, src) in &other.vertex_sources {
            self.vertex_sources.insert(*idx, *src);
        }
    }
}

/// Statistics about tracked history.
#[derive(Debug, Clone, Default)]
pub struct HistoryStatistics {
    pub modified_faces: usize,
    pub modified_edges: usize,
    pub modified_vertices: usize,
    pub generated_faces: usize,
    pub generated_edges: usize,
    pub generated_vertices: usize,
    pub deleted_faces: usize,
    pub deleted_edges: usize,
    pub deleted_vertices: usize,
}

impl HistoryStatistics {
    /// Total modified entities.
    pub fn total_modified(&self) -> usize {
        self.modified_faces + self.modified_edges + self.modified_vertices
    }

    /// Total generated entities.
    pub fn total_generated(&self) -> usize {
        self.generated_faces + self.generated_edges + self.generated_vertices
    }

    /// Total deleted entities.
    pub fn total_deleted(&self) -> usize {
        self.deleted_faces + self.deleted_edges + self.deleted_vertices
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BooleanHistory (backward compatible)
// ─────────────────────────────────────────────────────────────────────────────

/// Per-face origin map for a boolean operation result.
///
/// `face_origins[i]` gives the origin of `result_brep.solids[0].shells[0].faces[i]`.
#[derive(Debug, Clone)]
pub struct BooleanHistory {
    pub face_origins: Vec<FaceOrigin>,
    /// Per-edge origin map. `edge_origins[i]` gives the origin of `result_brep.edges[i]`.
    ///
    /// Empty when edge history was not requested (standard `boolean_op_with_history` path
    /// does not yet populate this; use `boolean_op_with_full_history` for edge tracking).
    pub edge_origins: Vec<EdgeOrigin>,
    /// Per-vertex origin map. `vertex_origins[i]` gives the origin of `result_brep.vertices[i]`.
    ///
    /// Empty when vertex history was not requested.
    pub vertex_origins: Vec<VertexOrigin>,
    /// Per-shell aggregate origin map. Flattened in the same order as `result_brep.solids[*].shells[*]`.
    pub shell_origins: Vec<ShellOrigin>,
    /// Per-solid aggregate origin map. `solid_origins[i]` gives the origin of `result_brep.solids[i]`.
    pub solid_origins: Vec<SolidOrigin>,
    /// Comprehensive history tracker for advanced queries.
    pub tracker: HistoryTracker,
    /// Deleted entities from input A (indices).
    pub deleted_from_a: Vec<usize>,
    /// Deleted entities from input B (indices).
    pub deleted_from_b: Vec<usize>,
    /// Deletion reasons, indexed by entity type and index.
    pub deletion_reasons: HashMap<(EntityType, usize), DeletionReason>,
}

/// Report produced when propagating persistent names through a boolean result.
///
/// This captures which source names could not be mapped into the result and
/// which names had to be duplicated because a single source entity generated
/// multiple result entities.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BooleanNamingPropagationReport {
    pub dropped_from_a: Vec<String>,
    pub dropped_from_b: Vec<String>,
    pub duplicated_from_a: Vec<String>,
    pub duplicated_from_b: Vec<String>,
}

impl BooleanHistory {
    /// Create a new empty BooleanHistory.
    pub fn new() -> Self {
        Self {
            face_origins: Vec::new(),
            edge_origins: Vec::new(),
            vertex_origins: Vec::new(),
            shell_origins: Vec::new(),
            solid_origins: Vec::new(),
            tracker: HistoryTracker::new(),
            deleted_from_a: Vec::new(),
            deleted_from_b: Vec::new(),
            deletion_reasons: HashMap::new(),
        }
    }

    /// Returns the origin of face `idx` in the result BRep.
    pub fn face_origin(&self, idx: usize) -> FaceOrigin {
        self.face_origins[idx]
    }

    /// Returns the origin of edge `idx` in the result BRep.
    /// Returns `None` if edge history was not recorded.
    pub fn edge_origin(&self, idx: usize) -> Option<EdgeOrigin> {
        self.edge_origins.get(idx).copied()
    }

    /// Returns the origin of vertex `idx` in the result BRep.
    /// Returns `None` if vertex history was not recorded.
    pub fn vertex_origin(&self, idx: usize) -> Option<VertexOrigin> {
        self.vertex_origins.get(idx).copied()
    }

    /// Returns the aggregate origin of shell `idx` in the flattened result BRep.
    pub fn shell_origin(&self, idx: usize) -> Option<ShellOrigin> {
        self.shell_origins.get(idx).copied()
    }

    /// Returns the aggregate origin of solid `idx` in the result BRep.
    pub fn solid_origin(&self, idx: usize) -> Option<SolidOrigin> {
        self.solid_origins.get(idx).copied()
    }

    /// Number of result faces tracked.
    pub fn len(&self) -> usize {
        self.face_origins.len()
    }

    pub fn is_empty(&self) -> bool {
        self.face_origins.is_empty()
    }

    /// How many result faces came from solid A.
    pub fn count_from_a(&self) -> usize {
        self.face_origins
            .iter()
            .filter(|o| matches!(o, FaceOrigin::FromA(_)))
            .count()
    }

    /// How many result faces came from solid B.
    pub fn count_from_b(&self) -> usize {
        self.face_origins
            .iter()
            .filter(|o| matches!(o, FaceOrigin::FromB(_)))
            .count()
    }

    /// How many result edges came from solid A (including splits).
    pub fn edge_count_from_a(&self) -> usize {
        self.edge_origins
            .iter()
            .filter(|o| matches!(o, EdgeOrigin::FromA(_) | EdgeOrigin::SplitFromA(_)))
            .count()
    }

    /// How many result edges came from solid B (including splits).
    pub fn edge_count_from_b(&self) -> usize {
        self.edge_origins
            .iter()
            .filter(|o| matches!(o, EdgeOrigin::FromB(_) | EdgeOrigin::SplitFromB(_)))
            .count()
    }

    /// How many result edges were generated at the intersection.
    pub fn edge_count_generated(&self) -> usize {
        self.edge_origins
            .iter()
            .filter(|o| matches!(o, EdgeOrigin::Generated))
            .count()
    }

    /// How many result solids contain contributions from both inputs and/or generated topology.
    pub fn solid_count_mixed(&self) -> usize {
        self.solid_origins
            .iter()
            .filter(|o| matches!(o, SolidOrigin::Mixed))
            .count()
    }

    // ── OCCT-style Query Methods ─────────────────────────────────────────────────

    /// Returns true if any entities were modified.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::HasModified()`.
    pub fn has_modified(&self) -> bool {
        self.tracker.has_modified()
    }

    /// Returns true if any entities were generated.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::HasGenerated()`.
    pub fn has_generated(&self) -> bool {
        self.tracker.has_generated()
    }

    /// Returns true if any entities were deleted.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::HasDeleted()`.
    pub fn has_deleted(&self) -> bool {
        self.tracker.has_deleted()
    }

    /// Returns true if a face from the input was deleted.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::IsDeleted()`.
    pub fn is_face_deleted(&self, source_face_idx: usize, from_a: bool) -> bool {
        if from_a {
            self.deleted_from_a.contains(&source_face_idx)
        } else {
            self.deleted_from_b.contains(&source_face_idx)
        }
    }

    /// Returns true if a face was modified during the operation.
    pub fn is_face_modified(&self, source_face_idx: usize, from_a: bool) -> bool {
        for (idx, origin) in self.face_origins.iter().enumerate() {
            match origin {
                FaceOrigin::FromA(src_idx) if *src_idx == source_face_idx && from_a => {
                    return self.tracker.is_face_modified(idx) || self.tracker.modified_faces(*src_idx).len() > 0;
                }
                FaceOrigin::FromB(src_idx) if *src_idx == source_face_idx && !from_a => {
                    return self.tracker.is_face_modified(idx) || self.tracker.modified_faces(*src_idx).len() > 0;
                }
                _ => {}
            }
        }
        false
    }

    /// Get the result faces that came from a source face.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::Modified()`.
    pub fn modified_faces(&self, source_idx: usize, from_a: bool) -> Vec<usize> {
        self.face_origins
            .iter()
            .enumerate()
            .filter_map(|(result_idx, origin)| match (origin, from_a) {
                (FaceOrigin::FromA(src), true) if *src == source_idx => Some(result_idx),
                (FaceOrigin::FromB(src), false) if *src == source_idx => Some(result_idx),
                _ => None,
            })
            .collect()
    }

    /// Get the result edges that came from a source edge (including splits).
    pub fn modified_edges(&self, source_idx: usize, from_a: bool) -> Vec<usize> {
        self.edge_origins
            .iter()
            .enumerate()
            .filter_map(|(result_idx, origin)| match (origin, from_a) {
                (EdgeOrigin::FromA(src), true) if *src == source_idx => Some(result_idx),
                (EdgeOrigin::FromB(src), false) if *src == source_idx => Some(result_idx),
                (EdgeOrigin::SplitFromA(src), true) if *src == source_idx => Some(result_idx),
                (EdgeOrigin::SplitFromB(src), false) if *src == source_idx => Some(result_idx),
                _ => None,
            })
            .collect()
    }

    /// Get the generated faces (intersection faces, etc.).
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::Generated()`.
    pub fn generated_faces(&self) -> Vec<usize> {
        self.face_origins
            .iter()
            .enumerate()
            .filter_map(|(idx, origin)| match origin {
                FaceOrigin::Generated => Some(idx),
                _ => None,
            })
            .collect()
    }

    /// Get the generated edges (intersection edges).
    pub fn generated_edges(&self) -> Vec<usize> {
        self.edge_origins
            .iter()
            .enumerate()
            .filter_map(|(idx, origin)| match origin {
                EdgeOrigin::Generated => Some(idx),
                _ => None,
            })
            .collect()
    }

    /// Get the generated vertices (intersection vertices).
    pub fn generated_vertices(&self) -> Vec<usize> {
        self.vertex_origins
            .iter()
            .enumerate()
            .filter_map(|(idx, origin)| match origin {
                VertexOrigin::Intersection => Some(idx),
                _ => None,
            })
            .collect()
    }

    /// Get the original entity for a result face.
    /// Returns (source_index, from_a) if the face came from an input.
    pub fn get_face_source(&self, result_idx: usize) -> Option<(usize, bool)> {
        match self.face_origins.get(result_idx)? {
            FaceOrigin::FromA(src) => Some((*src, true)),
            FaceOrigin::FromB(src) => Some((*src, false)),
            FaceOrigin::Generated => None,
        }
    }

    /// Get the original entity for a result edge.
    /// Returns (source_index, from_a) if the edge came from an input.
    pub fn get_edge_source(&self, result_idx: usize) -> Option<(usize, bool)> {
        match self.edge_origins.get(result_idx)? {
            EdgeOrigin::FromA(src) => Some((*src, true)),
            EdgeOrigin::FromB(src) => Some((*src, false)),
            EdgeOrigin::SplitFromA(src) => Some((*src, true)),
            EdgeOrigin::SplitFromB(src) => Some((*src, false)),
            EdgeOrigin::Generated => None,
        }
    }

    /// Get the original entity for a result vertex.
    /// Returns (source_index, from_a) if the vertex came from an input.
    pub fn get_vertex_source(&self, result_idx: usize) -> Option<(usize, bool)> {
        match self.vertex_origins.get(result_idx)? {
            VertexOrigin::FromA(src) => Some((*src, true)),
            VertexOrigin::FromB(src) => Some((*src, false)),
            VertexOrigin::Intersection => None,
        }
    }

    /// Populate the tracker from the origin vectors.
    /// Call this after setting up face_origins, edge_origins, vertex_origins.
    pub fn populate_tracker(&mut self) {
        // Populate face modifications and generations.
        for (result_idx, origin) in self.face_origins.iter().enumerate() {
            match origin {
                FaceOrigin::FromA(src_idx) => {
                    self.tracker.record_face_modified_multi(*src_idx, vec![result_idx], InputSource::A);
                }
                FaceOrigin::FromB(src_idx) => {
                    self.tracker.record_face_modified_multi(*src_idx, vec![result_idx], InputSource::B);
                }
                FaceOrigin::Generated => {
                    self.tracker.record_face_generated(result_idx, GenerationCause::Intersection);
                }
            }
        }

        // Populate edge modifications and generations.
        for (result_idx, origin) in self.edge_origins.iter().enumerate() {
            match origin {
                EdgeOrigin::FromA(src_idx) => {
                    self.tracker.record_edge_modified(*src_idx, vec![result_idx], InputSource::A, ModificationType::Preserved);
                }
                EdgeOrigin::FromB(src_idx) => {
                    self.tracker.record_edge_modified(*src_idx, vec![result_idx], InputSource::B, ModificationType::Preserved);
                }
                EdgeOrigin::SplitFromA(src_idx) => {
                    // Find all edges that split from this source.
                    let splits: Vec<usize> = self.edge_origins
                        .iter()
                        .enumerate()
                        .filter_map(|(idx, o)| match o {
                            EdgeOrigin::SplitFromA(s) if *s == *src_idx => Some(idx),
                            _ => None,
                        })
                        .collect();
                    self.tracker.record_edge_modified(*src_idx, splits, InputSource::A, ModificationType::Split);
                }
                EdgeOrigin::SplitFromB(src_idx) => {
                    let splits: Vec<usize> = self.edge_origins
                        .iter()
                        .enumerate()
                        .filter_map(|(idx, o)| match o {
                            EdgeOrigin::SplitFromB(s) if *s == *src_idx => Some(idx),
                            _ => None,
                        })
                        .collect();
                    self.tracker.record_edge_modified(*src_idx, splits, InputSource::B, ModificationType::Split);
                }
                EdgeOrigin::Generated => {
                    self.tracker.record_edge_generated(result_idx, GenerationCause::Intersection);
                }
            }
        }

        // Populate vertex modifications and generations.
        for (result_idx, origin) in self.vertex_origins.iter().enumerate() {
            match origin {
                VertexOrigin::FromA(src_idx) => {
                    self.tracker.record_vertex_modified(*src_idx, vec![result_idx], InputSource::A);
                }
                VertexOrigin::FromB(src_idx) => {
                    self.tracker.record_vertex_modified(*src_idx, vec![result_idx], InputSource::B);
                }
                VertexOrigin::Intersection => {
                    self.tracker.record_vertex_generated(result_idx, GenerationCause::Intersection);
                }
            }
        }

        // Populate deletions.
        for &deleted_idx in &self.deleted_from_a {
            let reason = self.deletion_reasons
                .get(&(EntityType::Face, deleted_idx))
                .cloned()
                .unwrap_or(DeletionReason::BooleanOperation);
            self.tracker.record_face_deleted_with_source(deleted_idx, reason, InputSource::A);
        }
        for &deleted_idx in &self.deleted_from_b {
            let reason = self.deletion_reasons
                .get(&(EntityType::Face, deleted_idx))
                .cloned()
                .unwrap_or(DeletionReason::BooleanOperation);
            self.tracker.record_face_deleted_with_source(deleted_idx, reason, InputSource::B);
        }
    }

    /// Convert this boolean history to naming events for integration with PersistentNamingEngine.
    ///
    /// This creates a sequence of naming events that represent the entity mappings
    /// captured by this history. The events can be applied to a PersistentNamingEngine
    /// to update its naming context.
    ///
    /// # Arguments
    /// * `_result_brep` - The result BRep (used for entity counts in future extensions).
    /// * `_entity_count_before_a` - Total entity count in solid A before the operation (for future use).
    /// * `_entity_count_before_b` - Total entity count in solid B before the operation (for future use).
    ///
    /// # Returns
    /// A vector of naming events representing the boolean operation's effect on naming.
    pub fn to_naming_events(
        &self,
        _result_brep: &BRep,
        _entity_count_before_a: usize,
        _entity_count_before_b: usize,
    ) -> Vec<NamingEvent> {
        let mut events = Vec::new();

        // Process face origins - create propagation events for faces from A/B.
        for (result_idx, origin) in self.face_origins.iter().enumerate() {
            let result_entity_id = Self::face_entity_id(result_idx);
            match origin {
                FaceOrigin::FromA(source_idx) => {
                    let source_entity_id = Self::face_entity_id(*source_idx);
                    // We'll create a generic propagation event.
                    // The persistent ID will be assigned by the engine.
                    events.push(NamingEvent::Propagated {
                        from_entity: source_entity_id,
                        to_entity: result_entity_id,
                        persistent_id: PersistentId::NULL, // Placeholder - engine will assign
                    });
                }
                FaceOrigin::FromB(source_idx) => {
                    let source_entity_id = Self::face_entity_id(*source_idx);
                    events.push(NamingEvent::Propagated {
                        from_entity: source_entity_id,
                        to_entity: result_entity_id,
                        persistent_id: PersistentId::NULL,
                    });
                }
                FaceOrigin::Generated => {
                    // Generated faces get new names.
                    events.push(NamingEvent::Assigned {
                        entity_id: result_entity_id,
                        persistent_id: PersistentId::NULL,
                    });
                }
            }
        }

        // Process edge origins - handle splits specially.
        let mut split_sources: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
        for (result_idx, origin) in self.edge_origins.iter().enumerate() {
            let result_entity_id = Self::edge_entity_id(result_idx);
            match origin {
                EdgeOrigin::FromA(source_idx) => {
                    let source_entity_id = Self::edge_entity_id(*source_idx);
                    events.push(NamingEvent::Propagated {
                        from_entity: source_entity_id,
                        to_entity: result_entity_id,
                        persistent_id: PersistentId::NULL,
                    });
                }
                EdgeOrigin::FromB(source_idx) => {
                    let source_entity_id = Self::edge_entity_id(*source_idx);
                    events.push(NamingEvent::Propagated {
                        from_entity: source_entity_id,
                        to_entity: result_entity_id,
                        persistent_id: PersistentId::NULL,
                    });
                }
                EdgeOrigin::SplitFromA(source_idx) => {
                    split_sources.entry(*source_idx).or_default().push(result_idx);
                }
                EdgeOrigin::SplitFromB(source_idx) => {
                    split_sources.entry(*source_idx).or_default().push(result_idx);
                }
                EdgeOrigin::Generated => {
                    events.push(NamingEvent::Assigned {
                        entity_id: result_entity_id,
                        persistent_id: PersistentId::NULL,
                    });
                }
            }
        }

        // Create split events for edges that were split.
        for (source_idx, target_indices) in split_sources {
            let source_entity_id = Self::edge_entity_id(source_idx);
            let target_entity_ids: Vec<u64> = target_indices.iter()
                .map(|&idx| Self::edge_entity_id(idx))
                .collect();
            // Create placeholder persistent IDs - engine will assign actual IDs.
            let target_persistent_ids: Vec<PersistentId> = target_indices.iter()
                .map(|_| PersistentId::NULL)
                .collect();
            events.push(NamingEvent::Split {
                source_entity: source_entity_id,
                target_entities: target_entity_ids,
                source_persistent_id: PersistentId::NULL,
                target_persistent_ids,
            });
        }

        // Process vertex origins.
        for (result_idx, origin) in self.vertex_origins.iter().enumerate() {
            let result_entity_id = Self::vertex_entity_id(result_idx);
            match origin {
                VertexOrigin::FromA(source_idx) => {
                    let source_entity_id = Self::vertex_entity_id(*source_idx);
                    events.push(NamingEvent::Propagated {
                        from_entity: source_entity_id,
                        to_entity: result_entity_id,
                        persistent_id: PersistentId::NULL,
                    });
                }
                VertexOrigin::FromB(source_idx) => {
                    let source_entity_id = Self::vertex_entity_id(*source_idx);
                    events.push(NamingEvent::Propagated {
                        from_entity: source_entity_id,
                        to_entity: result_entity_id,
                        persistent_id: PersistentId::NULL,
                    });
                }
                VertexOrigin::Intersection => {
                    events.push(NamingEvent::Assigned {
                        entity_id: result_entity_id,
                        persistent_id: PersistentId::NULL,
                    });
                }
            }
        }

        events
    }

    /// Helper to compute entity ID for a face index.
    /// Uses a simple encoding: entity_type * 1e9 + index.
    fn face_entity_id(idx: usize) -> u64 {
        1_000_000_000 + idx as u64
    }

    /// Helper to compute entity ID for an edge index.
    fn edge_entity_id(idx: usize) -> u64 {
        2_000_000_000 + idx as u64
    }

    /// Helper to compute entity ID for a vertex index.
    fn vertex_entity_id(idx: usize) -> u64 {
        3_000_000_000 + idx as u64
    }

    /// Apply this boolean history to a PersistentNamingEngine.
    ///
    /// This method begins an operation on the engine, applies all naming events,
    /// and finalizes the operation with statistics.
    ///
    /// # Arguments
    /// * `engine` - The naming engine to update.
    /// * `result_brep` - The result BRep.
    /// * `operation_type` - The type of boolean operation.
    /// * `entity_count_before_a` - Entity count in solid A before operation.
    /// * `entity_count_before_b` - Entity count in solid B before operation.
    ///
    /// # Returns
    /// The operation ID assigned by the engine.
    pub fn apply_to_naming_engine(
        &self,
        engine: &mut PersistentNamingEngine,
        result_brep: &BRep,
        operation_type: BooleanOperationType,
        entity_count_before_a: usize,
        entity_count_before_b: usize,
    ) -> rcad_kernel::persistent_naming::OperationId {
        let op_type = operation_type.into();

        let op_id = engine.begin_operation(op_type, None);

        let events = self.to_naming_events(
            result_brep,
            entity_count_before_a,
            entity_count_before_b,
        );

        // Apply events to the engine's context.
        for event in &events {
            engine.apply_and_track(event.clone());
        }

        let stats = OperationStats {
            entity_count_before: entity_count_before_a + entity_count_before_b,
            entity_count_after: result_brep.vertices.len() + result_brep.edges.len() +
                result_brep.solids.iter()
                    .flat_map(|s| s.shells.iter())
                    .map(|sh| sh.faces.len())
                    .sum::<usize>(),
            names_preserved: events.iter().filter(|e| matches!(e, NamingEvent::Propagated { .. })).count(),
            names_lost: 0, // Would need before context to determine
            names_generated: events.iter().filter(|e| matches!(e, NamingEvent::Assigned { .. })).count(),
            conflicts_resolved: 0,
        };

        engine.finalize_operation(stats);

        op_id
    }

    /// Propagate source persistent names through this boolean history into the
    /// `result_brep` topology.
    ///
    /// Face, edge, vertex, and baseline solid names from both inputs are mapped
    /// onto result entities according to `face_origins`, `edge_origins`,
    /// `vertex_origins`, and `solid_origins`.
    ///
    /// When one source entity produces multiple result entities (for example a
    /// split edge), the original name is bound to the first result entity and
    /// deterministic suffixed variants (`name@1`, `name@2`, ...) are bound to the
    /// remaining ones.
    pub fn propagate_persistent_naming(
        &self,
        result_brep: &BRep,
        names_a: &PersistentNamingHooks,
        names_b: &PersistentNamingHooks,
    ) -> (PersistentNamingHooks, BooleanNamingPropagationReport) {
        let mut out = PersistentNamingHooks::new();
        let mut report = BooleanNamingPropagationReport::default();

        self.propagate_from_source(
            &mut out,
            names_a,
            InputSide::A,
            &mut report.dropped_from_a,
            &mut report.duplicated_from_a,
        );
        self.propagate_from_source(
            &mut out,
            names_b,
            InputSide::B,
            &mut report.dropped_from_b,
            &mut report.duplicated_from_b,
        );

        out.retain_valid_for_brep(result_brep);
        (out, report)
    }

    fn propagate_from_source(
        &self,
        out: &mut PersistentNamingHooks,
        source: &PersistentNamingHooks,
        side: InputSide,
        dropped: &mut Vec<String>,
        duplicated: &mut Vec<String>,
    ) {
        for (name, target_ref) in source.iter() {
            let matches = self.matching_result_entities(side, target_ref);
            if matches.is_empty() {
                dropped.push(name.to_string());
                continue;
            }
            if matches.len() > 1 {
                duplicated.push(name.to_string());
            }
            for (suffix_idx, result_ref) in matches.into_iter().enumerate() {
                let bound_name = if suffix_idx == 0 {
                    name.to_string()
                } else {
                    format!("{name}@{suffix_idx}")
                };
                bind_unique(out, bound_name, result_ref);
            }
        }
    }

    fn matching_result_entities(&self, side: InputSide, target_ref: TopoEntityRef) -> Vec<TopoEntityRef> {
        match target_ref {
            TopoEntityRef::Face(source_idx) => self
                .face_origins
                .iter()
                .enumerate()
                .filter_map(|(result_idx, origin)| match (side, origin) {
                    (InputSide::A, FaceOrigin::FromA(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Face(result_idx))
                    }
                    (InputSide::B, FaceOrigin::FromB(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Face(result_idx))
                    }
                    _ => None,
                })
                .collect(),
            TopoEntityRef::Edge(source_idx) => self
                .edge_origins
                .iter()
                .enumerate()
                .filter_map(|(result_idx, origin)| match (side, origin) {
                    (InputSide::A, EdgeOrigin::FromA(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Edge(result_idx))
                    }
                    (InputSide::A, EdgeOrigin::SplitFromA(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Edge(result_idx))
                    }
                    (InputSide::B, EdgeOrigin::FromB(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Edge(result_idx))
                    }
                    (InputSide::B, EdgeOrigin::SplitFromB(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Edge(result_idx))
                    }
                    _ => None,
                })
                .collect(),
            TopoEntityRef::Vertex(source_idx) => self
                .vertex_origins
                .iter()
                .enumerate()
                .filter_map(|(result_idx, origin)| match (side, origin) {
                    (InputSide::A, VertexOrigin::FromA(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Vertex(result_idx))
                    }
                    (InputSide::B, VertexOrigin::FromB(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Vertex(result_idx))
                    }
                    _ => None,
                })
                .collect(),
            TopoEntityRef::Solid(source_idx) => {
                if source_idx != 0 {
                    return Vec::new();
                }
                self.solid_origins
                    .iter()
                    .enumerate()
                    .filter_map(|(result_idx, origin)| match (side, origin) {
                        (InputSide::A, SolidOrigin::FromA) => Some(TopoEntityRef::Solid(result_idx)),
                        (InputSide::B, SolidOrigin::FromB) => Some(TopoEntityRef::Solid(result_idx)),
                        _ => None,
                    })
                    .collect()
            }
        }
    }
}

impl Default for BooleanHistory {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputSide {
    A,
    B,
}

fn bind_unique(out: &mut PersistentNamingHooks, preferred_name: String, target_ref: TopoEntityRef) {
    if out.resolve(&preferred_name).is_none() {
        out.bind(preferred_name, target_ref);
        return;
    }
    let mut suffix_idx = 1usize;
    loop {
        let candidate = format!("{preferred_name}@{suffix_idx}");
        if out.resolve(&candidate).is_none() {
            out.bind(candidate, target_ref);
            return;
        }
        suffix_idx += 1;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Multi-Operation History Chain
// ─────────────────────────────────────────────────────────────────────────────

/// Chain of histories across multiple operations.
///
/// This enables ancestor/descendant queries across multiple boolean operations,
/// similar to OCCT's naming history tracking.
#[derive(Debug, Clone, Default)]
pub struct HistoryChain {
    /// Operations in the chain.
    operations: Vec<HistoryTracker>,
    /// Labels for each operation.
    labels: Vec<Option<String>>,
}

impl HistoryChain {
    /// Create a new empty history chain.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an operation to the chain.
    pub fn push(&mut self, history: HistoryTracker, label: Option<String>) {
        self.operations.push(history);
        self.labels.push(label);
    }

    /// Number of operations in the chain.
    pub fn len(&self) -> usize {
        self.operations.len()
    }

    pub fn is_empty(&self) -> bool {
        self.operations.is_empty()
    }

    /// Get a specific operation by index.
    pub fn get(&self, index: usize) -> Option<&HistoryTracker> {
        self.operations.get(index)
    }

    /// Get the label for a specific operation.
    pub fn label(&self, index: usize) -> Option<&str> {
        self.labels.get(index).and_then(|l| l.as_deref())
    }

    /// Trace the ancestry of an entity through all operations.
    ///
    /// Returns a vector of (operation_index, entity_index, entity_type) showing
    /// the lineage from the given entity back to its original source.
    pub fn trace_ancestry(&self, entity_type: EntityType, entity_index: usize) -> Vec<(usize, usize, EntityType)> {
        let mut result = Vec::new();
        let current_type = entity_type;
        let mut current_index = entity_index;

        // Walk backwards through operations.
        for (op_idx, tracker) in self.operations.iter().enumerate().rev() {
            result.push((op_idx, current_index, current_type));

            if let Some((source, source_idx)) = tracker.get_source(current_type, current_index) {
                if source == InputSource::Generated {
                    // Stop at generated entities - they have no further ancestry.
                    break;
                }
                current_index = source_idx;
            } else {
                // No source found, stop tracing.
                break;
            }
        }

        result
    }

    /// Find all descendants of an entity through all operations.
    ///
    /// Returns a vector of (operation_index, entity_indices) showing
    /// what entities were derived from the given source entity.
    pub fn trace_descendants(&self, entity_type: EntityType, source_index: usize) -> Vec<(usize, Vec<usize>)> {
        let mut result = Vec::new();
        let mut current_indices = vec![source_index];

        // Walk forwards through operations.
        for (op_idx, tracker) in self.operations.iter().enumerate() {
            let mut next_indices = Vec::new();

            for &idx in &current_indices {
                let descendants = tracker.get_results(entity_type, idx);
                next_indices.extend(descendants);
            }

            if !next_indices.is_empty() {
                result.push((op_idx, next_indices.clone()));
                current_indices = next_indices;
            }
        }

        result
    }

    /// Check if an entity is deleted in any operation.
    pub fn is_deleted_any(&self, entity_type: EntityType, entity_index: usize) -> bool {
        self.operations.iter().any(|tracker| tracker.is_deleted(entity_index, entity_type))
    }

    /// Get the deletion reason for an entity.
    pub fn get_deletion_reason(&self, entity_type: EntityType, entity_index: usize) -> Option<&DeletionReason> {
        for tracker in &self.operations {
            if let Some(record) = tracker.deletion_record(entity_index, entity_type) {
                return Some(&record.reason);
            }
        }
        None
    }

    /// Get comprehensive statistics across all operations.
    pub fn statistics(&self) -> ChainStatistics {
        let mut stats = ChainStatistics::default();

        for tracker in &self.operations {
            let op_stats = tracker.statistics();
            stats.total_modified_faces += op_stats.modified_faces;
            stats.total_modified_edges += op_stats.modified_edges;
            stats.total_modified_vertices += op_stats.modified_vertices;
            stats.total_generated_faces += op_stats.generated_faces;
            stats.total_generated_edges += op_stats.generated_edges;
            stats.total_generated_vertices += op_stats.generated_vertices;
            stats.total_deleted_faces += op_stats.deleted_faces;
            stats.total_deleted_edges += op_stats.deleted_edges;
            stats.total_deleted_vertices += op_stats.deleted_vertices;
        }

        stats.operation_count = self.operations.len();
        stats
    }

    /// Clear all operations.
    pub fn clear(&mut self) {
        self.operations.clear();
        self.labels.clear();
    }
}

/// Statistics across multiple operations.
#[derive(Debug, Clone, Default)]
pub struct ChainStatistics {
    pub operation_count: usize,
    pub total_modified_faces: usize,
    pub total_modified_edges: usize,
    pub total_modified_vertices: usize,
    pub total_generated_faces: usize,
    pub total_generated_edges: usize,
    pub total_generated_vertices: usize,
    pub total_deleted_faces: usize,
    pub total_deleted_edges: usize,
    pub total_deleted_vertices: usize,
}

impl ChainStatistics {
    pub fn total_modified(&self) -> usize {
        self.total_modified_faces + self.total_modified_edges + self.total_modified_vertices
    }

    pub fn total_generated(&self) -> usize {
        self.total_generated_faces + self.total_generated_edges + self.total_generated_vertices
    }

    pub fn total_deleted(&self) -> usize {
        self.total_deleted_faces + self.total_deleted_edges + self.total_deleted_vertices
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::{PrimitiveSolid, TopoEntityRef};

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    // ── HistoryTracker Tests ─────────────────────────────────────────────────────

    #[test]
    fn tracker_starts_empty() {
        let tracker = HistoryTracker::new();
        assert!(!tracker.has_modified());
        assert!(!tracker.has_generated());
        assert!(!tracker.has_deleted());
    }

    #[test]
    fn tracker_records_face_modification() {
        let mut tracker = HistoryTracker::new();
        tracker.record_face_modified(0, 1);

        assert!(tracker.has_modified());
        assert!(tracker.is_face_modified(0));
        assert_eq!(tracker.modified_faces(0), vec![1]);

        let source = tracker.get_source(EntityType::Face, 1);
        assert_eq!(source, Some((InputSource::A, 0)));
    }

    #[test]
    fn tracker_records_face_split() {
        let mut tracker = HistoryTracker::new();
        tracker.record_face_modified_multi(0, vec![1, 2, 3], InputSource::A);

        assert!(tracker.has_modified());
        let record = tracker.face_modification_record(0).unwrap();
        assert_eq!(record.modification_type, ModificationType::Split);
        assert_eq!(record.result_indices, vec![1, 2, 3]);
    }

    #[test]
    fn tracker_records_generated_entities() {
        let mut tracker = HistoryTracker::new();

        tracker.record_face_generated(10, GenerationCause::Intersection);
        tracker.record_edge_generated(20, GenerationCause::NewBoundary);
        tracker.record_vertex_generated(30, GenerationCause::Intersection);

        assert!(tracker.has_generated());
        assert!(tracker.is_face_generated(10));
        assert!(tracker.is_edge_generated(20));
        assert!(tracker.is_vertex_generated(30));

        assert!(!tracker.is_face_generated(11));
    }

    #[test]
    fn tracker_records_deletions() {
        let mut tracker = HistoryTracker::new();

        tracker.record_face_deleted(0, DeletionReason::BooleanOperation);
        tracker.record_edge_deleted(1, DeletionReason::Overlap);
        tracker.record_vertex_deleted(2, DeletionReason::Custom("Test".to_string()));

        assert!(tracker.has_deleted());
        assert!(tracker.is_face_deleted(0));
        assert!(tracker.is_edge_deleted(1));
        assert!(tracker.is_vertex_deleted(2));

        let record = tracker.deletion_record(0, EntityType::Face).unwrap();
        assert_eq!(record.reason, DeletionReason::BooleanOperation);

        let vertex_record = tracker.deletion_record(2, EntityType::Vertex).unwrap();
        if let DeletionReason::Custom(s) = &vertex_record.reason {
            assert_eq!(s, "Test");
        } else {
            panic!("Expected Custom deletion reason");
        }
    }

    #[test]
    fn tracker_count_by_source() {
        let mut tracker = HistoryTracker::new();

        tracker.record_face_modified_multi(0, vec![1, 2], InputSource::A);
        tracker.record_face_modified_multi(1, vec![3], InputSource::B);
        tracker.record_face_generated(4, GenerationCause::Intersection);

        assert_eq!(tracker.count_by_source(EntityType::Face, InputSource::A), 2);
        assert_eq!(tracker.count_by_source(EntityType::Face, InputSource::B), 1);
        assert_eq!(tracker.count_by_source(EntityType::Face, InputSource::Generated), 1);
    }

    #[test]
    fn tracker_merge() {
        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified(0, 1);

        let mut tracker2 = HistoryTracker::new();
        tracker2.record_face_modified(2, 3);
        tracker2.record_face_generated(4, GenerationCause::Intersection);

        tracker1.merge(&tracker2);

        assert!(tracker1.is_face_modified(0));
        assert!(tracker1.is_face_modified(2));
        assert!(tracker1.is_face_generated(4));
    }

    #[test]
    fn tracker_statistics() {
        let mut tracker = HistoryTracker::new();

        tracker.record_face_modified(0, 1);
        tracker.record_edge_modified(0, vec![1, 2], InputSource::A, ModificationType::Split);
        tracker.record_face_generated(10, GenerationCause::Intersection);
        tracker.record_face_deleted(20, DeletionReason::BooleanOperation);

        let stats = tracker.statistics();
        assert_eq!(stats.modified_faces, 1);
        assert_eq!(stats.modified_edges, 1);
        assert_eq!(stats.generated_faces, 1);
        assert_eq!(stats.deleted_faces, 1);
        assert_eq!(stats.total_modified(), 2);
        assert_eq!(stats.total_generated(), 1);
        assert_eq!(stats.total_deleted(), 1);
    }

    // ── BooleanHistory Tests ─────────────────────────────────────────────────────

    #[test]
    fn boolean_history_new() {
        let history = BooleanHistory::new();
        assert!(history.is_empty());
        assert!(!history.has_modified());
        assert!(!history.has_generated());
        assert!(!history.has_deleted());
    }

    #[test]
    fn boolean_history_populate_tracker() {
        let mut history = BooleanHistory::new();
        history.face_origins = vec![
            FaceOrigin::FromA(0),
            FaceOrigin::FromB(1),
            FaceOrigin::Generated,
        ];
        history.edge_origins = vec![
            EdgeOrigin::FromA(0),
            EdgeOrigin::SplitFromA(0),
            EdgeOrigin::Generated,
        ];
        history.vertex_origins = vec![
            VertexOrigin::FromA(0),
            VertexOrigin::Intersection,
        ];

        history.populate_tracker();

        assert!(history.has_modified());
        assert!(history.has_generated());

        // Check face queries.
        assert_eq!(history.modified_faces(0, true), vec![0]);
        assert_eq!(history.modified_faces(1, false), vec![1]);
        assert_eq!(history.generated_faces(), vec![2]);

        // Check edge queries.
        let modified_edges = history.modified_edges(0, true);
        assert_eq!(modified_edges.len(), 2);
        assert!(modified_edges.contains(&0));

        // Check vertex queries.
        assert_eq!(history.generated_vertices(), vec![1]);
    }

    #[test]
    fn boolean_history_get_source() {
        let mut history = BooleanHistory::new();
        history.face_origins = vec![
            FaceOrigin::FromA(0),
            FaceOrigin::FromB(5),
            FaceOrigin::Generated,
        ];

        assert_eq!(history.get_face_source(0), Some((0, true)));
        assert_eq!(history.get_face_source(1), Some((5, false)));
        assert_eq!(history.get_face_source(2), None);
    }

    #[test]
    fn boolean_history_deletions() {
        let mut history = BooleanHistory::new();
        history.deleted_from_a = vec![0, 1];
        history.deleted_from_b = vec![2];
        history.deletion_reasons.insert(
            (EntityType::Face, 0),
            DeletionReason::BooleanOperation,
        );

        assert!(history.is_face_deleted(0, true));
        assert!(history.is_face_deleted(1, true));
        assert!(history.is_face_deleted(2, false));
        assert!(!history.is_face_deleted(3, true));
    }

    #[test]
    fn propagate_persistent_naming_maps_face_edge_vertex_and_solid_origins() {
        let result_brep = unit_box();
        let mut history = BooleanHistory {
            face_origins: vec![FaceOrigin::FromA(0), FaceOrigin::FromB(1)],
            edge_origins: vec![EdgeOrigin::FromA(0), EdgeOrigin::SplitFromA(0), EdgeOrigin::FromB(1)],
            vertex_origins: vec![VertexOrigin::FromA(0), VertexOrigin::Intersection, VertexOrigin::FromB(1)],
            shell_origins: vec![],
            solid_origins: vec![SolidOrigin::FromA],
            tracker: HistoryTracker::new(),
            deleted_from_a: vec![],
            deleted_from_b: vec![],
            deletion_reasons: HashMap::new(),
        };

        let mut names_a = PersistentNamingHooks::new();
        names_a.bind("face_a", TopoEntityRef::Face(0));
        names_a.bind("edge_a", TopoEntityRef::Edge(0));
        names_a.bind("vertex_a", TopoEntityRef::Vertex(0));
        names_a.bind("solid_a", TopoEntityRef::Solid(0));

        let mut names_b = PersistentNamingHooks::new();
        names_b.bind("face_b", TopoEntityRef::Face(1));
        names_b.bind("edge_b", TopoEntityRef::Edge(1));
        names_b.bind("vertex_b", TopoEntityRef::Vertex(1));
        names_b.bind("solid_b", TopoEntityRef::Solid(0));

        let (result_names, report) = history.propagate_persistent_naming(&result_brep, &names_a, &names_b);

        assert_eq!(result_names.resolve("face_a"), Some(TopoEntityRef::Face(0)));
        assert_eq!(result_names.resolve("face_b"), Some(TopoEntityRef::Face(1)));
        assert_eq!(result_names.resolve("edge_a"), Some(TopoEntityRef::Edge(0)));
        assert_eq!(result_names.resolve("edge_a@1"), Some(TopoEntityRef::Edge(1)));
        assert_eq!(result_names.resolve("edge_b"), Some(TopoEntityRef::Edge(2)));
        assert_eq!(result_names.resolve("vertex_a"), Some(TopoEntityRef::Vertex(0)));
        assert_eq!(result_names.resolve("vertex_b"), Some(TopoEntityRef::Vertex(2)));
        assert_eq!(result_names.resolve("solid_a"), Some(TopoEntityRef::Solid(0)));
        assert_eq!(result_names.resolve("solid_b"), None);

        assert!(report.dropped_from_a.is_empty());
        assert_eq!(report.dropped_from_b, vec!["solid_b".to_string()]);
        assert_eq!(report.duplicated_from_a, vec!["edge_a".to_string()]);
        assert!(report.duplicated_from_b.is_empty());
    }

    #[test]
    fn propagate_persistent_naming_disambiguates_cross_input_name_collisions() {
        let result_brep = unit_box();
        let mut history = BooleanHistory {
            face_origins: vec![FaceOrigin::FromA(0), FaceOrigin::FromB(0)],
            edge_origins: vec![],
            vertex_origins: vec![],
            shell_origins: vec![],
            solid_origins: vec![],
            tracker: HistoryTracker::new(),
            deleted_from_a: vec![],
            deleted_from_b: vec![],
            deletion_reasons: HashMap::new(),
        };

        let mut names_a = PersistentNamingHooks::new();
        names_a.bind("shared_face", TopoEntityRef::Face(0));
        let mut names_b = PersistentNamingHooks::new();
        names_b.bind("shared_face", TopoEntityRef::Face(0));

        let (result_names, report) = history.propagate_persistent_naming(&result_brep, &names_a, &names_b);

        assert_eq!(result_names.resolve("shared_face"), Some(TopoEntityRef::Face(0)));
        assert_eq!(result_names.resolve("shared_face@1"), Some(TopoEntityRef::Face(1)));
        assert!(report.dropped_from_a.is_empty());
        assert!(report.dropped_from_b.is_empty());
    }

    // ── HistoryChain Tests ───────────────────────────────────────────────────────

    #[test]
    fn chain_starts_empty() {
        let chain = HistoryChain::new();
        assert!(chain.is_empty());
        assert_eq!(chain.len(), 0);
    }

    #[test]
    fn chain_push_and_get() {
        let mut chain = HistoryChain::new();

        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified(0, 1);
        chain.push(tracker1, Some("op1".to_string()));

        let mut tracker2 = HistoryTracker::new();
        tracker2.record_face_modified(1, 2);
        chain.push(tracker2, Some("op2".to_string()));

        assert_eq!(chain.len(), 2);
        assert_eq!(chain.label(0), Some("op1"));
        assert_eq!(chain.label(1), Some("op2"));
        assert!(chain.get(0).unwrap().is_face_modified(0));
        assert!(chain.get(1).unwrap().is_face_modified(1));
    }

    #[test]
    fn chain_trace_ancestry() {
        let mut chain = HistoryChain::new();

        // First operation: face 0 -> face 1.
        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified_multi(0, vec![1], InputSource::A);
        chain.push(tracker1, None);

        // Second operation: face 1 -> face 2.
        let mut tracker2 = HistoryTracker::new();
        tracker2.record_face_modified_multi(1, vec![2], InputSource::A);
        chain.push(tracker2, None);

        // Trace ancestry of face 2.
        let ancestry = chain.trace_ancestry(EntityType::Face, 2);
        assert_eq!(ancestry.len(), 2);

        // Should trace back through both operations.
        assert_eq!(ancestry[0].0, 1); // Second operation.
        assert_eq!(ancestry[0].1, 2); // Face 2.

        assert_eq!(ancestry[1].0, 0); // First operation.
        assert_eq!(ancestry[1].1, 1); // Face 1.
    }

    #[test]
    fn chain_trace_descendants() {
        let mut chain = HistoryChain::new();

        // First operation: face 0 -> faces 1, 2 (split).
        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified_multi(0, vec![1, 2], InputSource::A);
        chain.push(tracker1, None);

        // Second operation: faces 1, 2 -> faces 3, 4, 5.
        let mut tracker2 = HistoryTracker::new();
        tracker2.record_face_modified_multi(1, vec![3], InputSource::A);
        tracker2.record_face_modified_multi(2, vec![4, 5], InputSource::A);
        chain.push(tracker2, None);

        // Trace descendants of face 0.
        let descendants = chain.trace_descendants(EntityType::Face, 0);
        assert_eq!(descendants.len(), 2);

        // First operation: 0 -> [1, 2].
        assert_eq!(descendants[0].0, 0);
        assert_eq!(descendants[0].1.len(), 2);

        // Second operation: 1, 2 -> [3, 4, 5].
        assert_eq!(descendants[1].0, 1);
        assert_eq!(descendants[1].1.len(), 3);
    }

    #[test]
    fn chain_is_deleted_any() {
        let mut chain = HistoryChain::new();

        // First operation: delete face 0.
        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_deleted(0, DeletionReason::BooleanOperation);
        chain.push(tracker1, None);

        // Second operation: no deletions.
        let tracker2 = HistoryTracker::new();
        chain.push(tracker2, None);

        assert!(chain.is_deleted_any(EntityType::Face, 0));
        assert!(!chain.is_deleted_any(EntityType::Face, 1));
    }

    #[test]
    fn chain_statistics() {
        let mut chain = HistoryChain::new();

        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified(0, 1);
        tracker1.record_face_generated(2, GenerationCause::Intersection);
        chain.push(tracker1, None);

        let mut tracker2 = HistoryTracker::new();
        tracker2.record_edge_modified(0, vec![1, 2], InputSource::A, ModificationType::Split);
        tracker2.record_face_deleted(3, DeletionReason::BooleanOperation);
        chain.push(tracker2, None);

        let stats = chain.statistics();
        assert_eq!(stats.operation_count, 2);
        assert_eq!(stats.total_modified_faces, 1);
        assert_eq!(stats.total_modified_edges, 1);
        assert_eq!(stats.total_generated_faces, 1);
        assert_eq!(stats.total_deleted_faces, 1);
    }

    // ── DeletionReason Tests ─────────────────────────────────────────────────────

    #[test]
    fn deletion_reason_description() {
        assert_eq!(DeletionReason::BooleanOperation.description(), "Removed by boolean operation");
        assert_eq!(DeletionReason::OutsideResult.description(), "Outside result volume");
        assert_eq!(DeletionReason::Overlap.description(), "Overlapping geometry");
        assert_eq!(DeletionReason::Tolerance.description(), "Tolerance issues");
        assert_eq!(DeletionReason::Healing.description(), "Removed during healing");
        assert_eq!(
            DeletionReason::Custom("Custom reason".to_string()).description(),
            "Custom reason"
        );
    }

    // ── GenerationCause Tests ────────────────────────────────────────────────────

    #[test]
    fn generation_record_with_parents() {
        let mut tracker = HistoryTracker::new();
        tracker.record_edge_generated_with_parents(
            10,
            GenerationCause::Intersection,
            vec![0, 1], // Parent edges
        );

        let record = tracker.edge_generation_record(10).unwrap();
        assert_eq!(record.entity_index, 10);
        assert_eq!(record.cause, GenerationCause::Intersection);
        assert_eq!(record.parent_indices, vec![0, 1]);
    }
}

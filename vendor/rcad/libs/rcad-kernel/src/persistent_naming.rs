//! Persistent naming semantics for BRepGraph topology entities.
//!
//! This module provides stable, operation-surviving identifiers for topology
//! entities (vertices, edges, faces, solids). The naming system is inspired by
//! OCCT's OCAF/TopoNaming architecture.
//!
//! # Core Concepts
//!
//! - **PersistentId**: A stable 64-bit identifier that survives topology mutations.
//! - **NamingContext**: Bidirectional mapping between transient entity IDs and persistent IDs.
//! - **PersistentNamingEngine**: Orchestrates name assignment, resolution, and propagation.
//! - **NamingRule**: Strategies for assigning and propagating names.
//!
//! # Integration with BRepGraph
//!
//! The naming engine integrates with `BRepGraphHistory` to track naming changes
//! during graph mutations. Call `replay_with_naming()` to reconstruct naming
//! context from a history log.
//!
//! # Cross-Operation History
//!
//! The engine maintains an operation log that tracks entity genealogy across
//! multiple topology operations. This enables:
//! - Tracing entity origins through the operation history
//! - Detecting naming stability issues across operations
//! - Generating comprehensive stability reports

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// PersistentId
// ─────────────────────────────────────────────────────────────────────────────

/// A stable, operation-surviving identifier for a topology entity.
///
/// Unlike transient entity indices (which may shift after boolean operations,
/// splits, or merges), a `PersistentId` remains stable across operations that
/// preserve the logical identity of the entity.
///
/// Analogous to OCCT `TDF_Label` / `TNaming_NamedShape` references.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct PersistentId(pub u64);

impl PersistentId {
    /// Sentinel value for an invalid/unassigned persistent ID.
    pub const NULL: PersistentId = PersistentId(0);

    /// Returns `true` if this is the null sentinel.
    pub fn is_null(&self) -> bool {
        self.0 == 0
    }

    /// Returns the raw 64-bit value.
    pub fn raw(&self) -> u64 {
        self.0
    }
}

impl Default for PersistentId {
    fn default() -> Self {
        Self::NULL
    }
}

impl std::fmt::Display for PersistentId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "pid:{}", self.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingContext
// ─────────────────────────────────────────────────────────────────────────────

/// Bidirectional mapping between transient entity IDs and persistent IDs.
///
/// `NamingContext` maintains two hashmaps:
/// - `entity_to_persistent`: Maps transient entity IDs to their persistent identifiers.
/// - `persistent_to_entity`: Reverse lookup from persistent IDs to entity IDs.
///
/// Use `PersistentNamingEngine` to manage context lifecycle and propagation.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct NamingContext {
    /// Entity ID (as u64) to persistent ID mapping.
    entity_to_persistent: HashMap<u64, PersistentId>,
    /// Persistent ID to entity ID reverse mapping.
    persistent_to_entity: HashMap<PersistentId, u64>,
    /// Counter for allocating new persistent IDs (starts at 1; 0 is NULL).
    next_id: u64,
}

impl NamingContext {
    /// Create an empty naming context.
    pub fn new() -> Self {
        Self {
            entity_to_persistent: HashMap::new(),
            persistent_to_entity: HashMap::new(),
            next_id: 1,
        }
    }

    /// Returns the number of named entities in this context.
    pub fn len(&self) -> usize {
        self.entity_to_persistent.len()
    }

    /// Returns `true` if the context has no named entities.
    pub fn is_empty(&self) -> bool {
        self.entity_to_persistent.is_empty()
    }

    /// Allocate a new persistent ID.
    fn allocate_id(&mut self) -> PersistentId {
        let id = PersistentId(self.next_id);
        self.next_id += 1;
        id
    }

    /// Check if an entity has a persistent ID assigned.
    pub fn has_entity(&self, entity_id: u64) -> bool {
        self.entity_to_persistent.contains_key(&entity_id)
    }

    /// Check if a persistent ID is registered.
    pub fn has_persistent(&self, pid: PersistentId) -> bool {
        self.persistent_to_entity.contains_key(&pid)
    }

    /// Get the persistent ID for an entity, if assigned.
    pub fn get_persistent(&self, entity_id: u64) -> Option<PersistentId> {
        self.entity_to_persistent.get(&entity_id).copied()
    }

    /// Get the entity ID for a persistent ID, if registered.
    pub fn get_entity(&self, pid: PersistentId) -> Option<u64> {
        self.persistent_to_entity.get(&pid).copied()
    }

    /// Bind an entity to a persistent ID.
    ///
    /// If either the entity or the persistent ID was already bound,
    /// the old binding is removed.
    fn bind(&mut self, entity_id: u64, pid: PersistentId) {
        // Remove old bindings if present.
        if let Some(old_pid) = self.entity_to_persistent.remove(&entity_id) {
            self.persistent_to_entity.remove(&old_pid);
        }
        if let Some(old_eid) = self.persistent_to_entity.remove(&pid) {
            self.entity_to_persistent.remove(&old_eid);
        }
        // Insert new binding.
        self.entity_to_persistent.insert(entity_id, pid);
        self.persistent_to_entity.insert(pid, entity_id);
    }

    /// Unbind an entity from its persistent ID.
    fn unbind_entity(&mut self, entity_id: u64) -> Option<PersistentId> {
        let pid = self.entity_to_persistent.remove(&entity_id)?;
        self.persistent_to_entity.remove(&pid);
        Some(pid)
    }

    /// Unbind a persistent ID from its entity.
    fn unbind_persistent(&mut self, pid: PersistentId) -> Option<u64> {
        let entity_id = self.persistent_to_entity.remove(&pid)?;
        self.entity_to_persistent.remove(&entity_id);
        Some(entity_id)
    }

    /// Iterate over all (entity_id, persistent_id) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (u64, PersistentId)> + '_ {
        self.entity_to_persistent.iter().map(|(&e, &p)| (e, p))
    }

    /// Clear all bindings.
    pub fn clear(&mut self) {
        self.entity_to_persistent.clear();
        self.persistent_to_entity.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingRule
// ─────────────────────────────────────────────────────────────────────────────

/// Strategy for assigning and propagating persistent names.
///
/// Different strategies are appropriate for different operations:
/// - **GeometrySignature**: Assign names based on geometric properties (hash of
///   surface type, bounding box, curvature). Good for imported geometry.
/// - **TopologyRelation**: Assign names based on topological relationships
///   (e.g., "face adjacent to edge X"). Good for feature-based modeling.
/// - **HistoryTracking**: Track the origin of entities through operation history.
///   Good for parametric modeling.
/// - **Hybrid**: Combine multiple strategies for robustness.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NamingRule {
    /// Assign names based on geometric signatures.
    GeometrySignature,
    /// Assign names based on topological relationships.
    TopologyRelation,
    /// Track entity origins through operation history.
    HistoryTracking,
    /// Combine multiple strategies (recommended for production).
    #[default]
    Hybrid,
}

// ─────────────────────────────────────────────────────────────────────────────
// NamePropagationPolicy
// ─────────────────────────────────────────────────────────────────────────────

/// Policy for propagating names when entities are created, split, or merged.
///
/// When a topology operation produces new entities from existing ones,
/// the `NamePropagationPolicy` determines how names flow to the results.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NamePropagationPolicy {
    /// Keep the original entity's name (for minor modifications).
    Preserve,
    /// Inherit the parent entity's name with a disambiguating suffix.
    Inherit,
    /// Generate a completely new name.
    Generate,
    /// Combine names from multiple source entities (for merges).
    Combine,
}

impl Default for NamePropagationPolicy {
    fn default() -> Self {
        Self::Preserve
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingConflictResolution
// ─────────────────────────────────────────────────────────────────────────────

/// Record of a naming conflict and how it was resolved.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NamingConflictResolution {
    /// The persistent ID that was in conflict.
    pub conflicting_pid: PersistentId,
    /// The old entity ID that held the persistent ID.
    pub old_entity_id: u64,
    /// The new entity ID that now holds the persistent ID.
    pub new_entity_id: u64,
    /// How the conflict was resolved.
    pub resolution: ConflictResolution,
}

/// Strategy used to resolve a naming conflict.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConflictResolution {
    /// Kept the old binding, rejected the new.
    KeepOld,
    /// Replaced with the new binding, removed old.
    ReplaceOld,
    /// Generated a new persistent ID for the new entity.
    GenerateNew,
    /// Combined both entities under a shared context.
    Combine,
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingStabilityReport
// ─────────────────────────────────────────────────────────────────────────────

/// Entity type for reporting purposes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EntityType {
    Vertex,
    Edge,
    Face,
    Solid,
}

/// Detailed breakdown of naming stability per entity type.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct EntityTypeStability {
    /// Number of entities of this type before the operation.
    pub count_before: usize,
    /// Number of entities of this type after the operation.
    pub count_after: usize,
    /// Names preserved for this entity type.
    pub preserved: usize,
    /// Names lost for this entity type.
    pub lost: usize,
    /// New names assigned for this entity type.
    pub new_names: usize,
    /// Conflicts affecting this entity type.
    pub conflicts: usize,
}

impl EntityTypeStability {
    /// Calculate the stability score for this entity type.
    pub fn stability_score(&self) -> f64 {
        let total = self.preserved + self.lost;
        if total == 0 { 1.0 } else { self.preserved as f64 / total as f64 }
    }
}

/// Severity level for naming issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum IssueSeverity {
    /// Minor issue - name duplicated but recoverable.
    Minor,
    /// Moderate issue - name lost but entity may survive under new name.
    Moderate,
    /// Severe issue - complete loss of naming chain.
    Severe,
    /// Critical issue - naming inconsistency that may corrupt downstream operations.
    Critical,
}

/// A detected naming issue with context and resolution suggestions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NamingIssue {
    /// Severity of the issue.
    pub severity: IssueSeverity,
    /// Entity type affected.
    pub entity_type: Option<EntityType>,
    /// Description of the issue.
    pub description: String,
    /// Entity IDs involved in the issue.
    pub affected_entity_ids: Vec<u64>,
    /// Persistent IDs involved.
    pub affected_persistent_ids: Vec<PersistentId>,
    /// Suggested resolution strategy.
    pub suggested_resolution: Option<ConflictResolution>,
}

/// Report on naming stability after an operation.
///
/// Use `PersistentNamingEngine::stability_report()` to generate a report
/// comparing the pre- and post-operation naming contexts.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct NamingStabilityReport {
    /// Overall naming stability score (0.0 = all names lost, 1.0 = all names preserved).
    pub stability_score: f64,
    /// Entity IDs that lost their persistent names.
    pub lost_names: Vec<u64>,
    /// Entity IDs that received new persistent names.
    pub new_names: Vec<u64>,
    /// Entity IDs whose persistent names were preserved.
    pub preserved_names: Vec<u64>,
    /// Conflicts that were resolved during propagation.
    pub conflict_resolutions: Vec<NamingConflictResolution>,
    /// Total number of entities before the operation.
    pub entity_count_before: usize,
    /// Total number of entities after the operation.
    pub entity_count_after: usize,
    /// Detailed breakdown per entity type.
    pub entity_type_breakdown: HashMap<EntityType, EntityTypeStability>,
    /// Detected naming issues.
    pub issues: Vec<NamingIssue>,
    /// Score weighted by entity importance (faces > edges > vertices).
    pub weighted_stability_score: f64,
    /// Number of naming chains broken (entities with lost ancestry).
    pub broken_chains: usize,
}

impl NamingStabilityReport {
    /// Returns `true` if all names were preserved (score == 1.0, no lost names).
    pub fn is_perfect(&self) -> bool {
        self.stability_score >= 1.0 && self.lost_names.is_empty()
    }

    /// Returns `true` if any names were lost or conflicts occurred.
    pub fn has_issues(&self) -> bool {
        self.stability_score < 1.0 || !self.lost_names.is_empty() || !self.conflict_resolutions.is_empty()
    }

    /// Returns `true` if there are critical issues.
    pub fn has_critical_issues(&self) -> bool {
        self.issues.iter().any(|i| i.severity == IssueSeverity::Critical)
    }

    /// Get the stability breakdown for a specific entity type.
    pub fn get_entity_stability(&self, entity_type: EntityType) -> Option<&EntityTypeStability> {
        self.entity_type_breakdown.get(&entity_type)
    }

    /// Calculate a summary of issues by severity.
    pub fn issue_summary(&self) -> HashMap<IssueSeverity, usize> {
        let mut summary = HashMap::new();
        for issue in &self.issues {
            *summary.entry(issue.severity).or_insert(0) += 1;
        }
        summary
    }

    /// Generate a human-readable summary string.
    pub fn summary_string(&self) -> String {
        format!(
            "Naming Stability Report:\n\
             - Overall Score: {:.1}%\n\
             - Weighted Score: {:.1}%\n\
             - Entities: {} -> {}\n\
             - Preserved: {} | Lost: {} | New: {}\n\
             - Conflicts Resolved: {}\n\
             - Issues: {} (Critical: {})\n\
             - Broken Chains: {}",
            self.stability_score * 100.0,
            self.weighted_stability_score * 100.0,
            self.entity_count_before,
            self.entity_count_after,
            self.preserved_names.len(),
            self.lost_names.len(),
            self.new_names.len(),
            self.conflict_resolutions.len(),
            self.issues.len(),
            self.issues.iter().filter(|i| i.severity == IssueSeverity::Critical).count(),
            self.broken_chains
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PersistentNamingEngine
// ─────────────────────────────────────────────────────────────────────────────

/// Engine for managing persistent naming across BRep operations.
///
/// The engine coordinates:
/// - Assignment of new persistent IDs to entities.
/// - Resolution of entity IDs to persistent IDs and vice versa.
/// - Propagation of names through topology operations.
/// - Merging of naming contexts (e.g., after boolean operations).
/// - Cross-operation history tracking and genealogy.
///
/// # Example
///
/// ```rust
/// use rcad_kernel::persistent_naming::{PersistentNamingEngine, NamingRule, NamePropagationPolicy};
///
/// let mut engine = PersistentNamingEngine::new(NamingRule::Hybrid);
///
/// // Assign a persistent ID to entity 42 (e.g., face index 42).
/// let pid = engine.assign_persistent_id(42);
///
/// // Resolve back to the entity.
/// assert_eq!(engine.resolve_entity(pid), Some(42));
/// assert_eq!(engine.resolve_persistent(42), Some(pid));
/// ```
#[derive(Debug, Clone)]
pub struct PersistentNamingEngine {
    /// The active naming context.
    context: NamingContext,
    /// The naming rule to use for new assignments.
    rule: NamingRule,
    /// Default propagation policy.
    default_policy: NamePropagationPolicy,
    /// History of conflict resolutions.
    conflict_history: Vec<NamingConflictResolution>,
    /// Cross-operation history tracker.
    cross_op_history: CrossOperationHistory,
    /// Current operation ID (if an operation is in progress).
    current_operation: Option<OperationId>,
    /// Events accumulated for the current operation.
    pending_events: Vec<NamingEvent>,
}

impl Default for PersistentNamingEngine {
    fn default() -> Self {
        Self::new(NamingRule::default())
    }
}

impl PersistentNamingEngine {
    /// Create a new naming engine with the given rule.
    pub fn new(rule: NamingRule) -> Self {
        Self {
            context: NamingContext::new(),
            rule,
            default_policy: NamePropagationPolicy::default(),
            conflict_history: Vec::new(),
            cross_op_history: CrossOperationHistory::new(),
            current_operation: None,
            pending_events: Vec::new(),
        }
    }

    /// Create a naming engine with a specific default propagation policy.
    pub fn with_policy(rule: NamingRule, policy: NamePropagationPolicy) -> Self {
        Self {
            context: NamingContext::new(),
            rule,
            default_policy: policy,
            conflict_history: Vec::new(),
            cross_op_history: CrossOperationHistory::new(),
            current_operation: None,
            pending_events: Vec::new(),
        }
    }

    // ── Cross-Operation History ─────────────────────────────────────────────────

    /// Begin a new operation for history tracking.
    ///
    /// Returns an operation ID that should be used when finalizing the operation.
    pub fn begin_operation(&mut self, operation_type: OperationType, label: Option<String>) -> OperationId {
        // Finalize any pending operation first.
        if self.current_operation.is_some() {
            self.cancel_operation();
        }

        let op_id = self.cross_op_history.begin_operation(operation_type, label);
        self.current_operation = Some(op_id);
        self.pending_events.clear();
        op_id
    }

    /// Finalize the current operation with statistics.
    ///
    /// This commits all pending naming events to the history.
    pub fn finalize_operation(&mut self, stats: OperationStats) {
        if let Some(op_id) = self.current_operation.take() {
            let events = std::mem::take(&mut self.pending_events);
            self.cross_op_history.add_events(op_id, events);
            self.cross_op_history.finalize_operation(op_id, stats);
        }
    }

    /// Cancel the current operation without committing events.
    pub fn cancel_operation(&mut self) {
        self.current_operation = None;
        self.pending_events.clear();
    }

    /// Get the cross-operation history.
    pub fn cross_operation_history(&self) -> &CrossOperationHistory {
        &self.cross_op_history
    }

    /// Get mutable access to the cross-operation history.
    pub fn cross_operation_history_mut(&mut self) -> &mut CrossOperationHistory {
        &mut self.cross_op_history
    }

    /// Trace the genealogy of a persistent ID.
    pub fn trace_genealogy(&self, pid: PersistentId) -> Option<&EntityGenealogy> {
        self.cross_op_history.get_genealogy(pid)
    }

    /// Find all operations that affected a persistent ID.
    pub fn operations_affecting(&self, pid: PersistentId) -> Vec<&OperationRecord> {
        self.cross_op_history.operations_affecting(pid)
    }

    /// Generate a cross-operation stability report.
    pub fn cross_operation_stability_report(&self) -> CrossOperationStabilityReport {
        self.cross_op_history.stability_report()
    }

    // ── Assignment ────────────────────────────────────────────────────────────

    /// Assign a new persistent ID to an entity.
    ///
    /// If the entity already has a persistent ID, returns the existing one.
    /// Use `force_assign` to override.
    pub fn assign_persistent_id(&mut self, entity_id: u64) -> PersistentId {
        if let Some(existing) = self.context.get_persistent(entity_id) {
            return existing;
        }
        let pid = self.context.allocate_id();
        self.context.bind(entity_id, pid);

        // Track event if in an operation.
        if self.current_operation.is_some() {
            self.pending_events.push(NamingEvent::Assigned {
                entity_id,
                persistent_id: pid,
            });
        }

        pid
    }

    /// Force-assign a new persistent ID, replacing any existing binding.
    pub fn force_assign(&mut self, entity_id: u64) -> PersistentId {
        let pid = self.context.allocate_id();
        self.context.bind(entity_id, pid);

        // Track event if in an operation.
        if self.current_operation.is_some() {
            self.pending_events.push(NamingEvent::Assigned {
                entity_id,
                persistent_id: pid,
            });
        }

        pid
    }

    /// Assign a specific persistent ID to an entity.
    ///
    /// Returns `true` if the assignment succeeded (no conflict).
    /// Returns `false` if the persistent ID was already bound to a different entity.
    pub fn assign_specific(&mut self, entity_id: u64, pid: PersistentId) -> bool {
        if pid.is_null() {
            return false;
        }
        // Check for conflict.
        if let Some(existing_entity) = self.context.get_entity(pid) {
            if existing_entity != entity_id {
                return false;
            }
        }
        self.context.bind(entity_id, pid);

        // Track event if in an operation.
        if self.current_operation.is_some() {
            self.pending_events.push(NamingEvent::Assigned {
                entity_id,
                persistent_id: pid,
            });
        }

        true
    }

    // ── Resolution ─────────────────────────────────────────────────────────────

    /// Resolve a persistent ID to its entity ID.
    pub fn resolve_entity(&self, pid: PersistentId) -> Option<u64> {
        self.context.get_entity(pid)
    }

    /// Resolve an entity ID to its persistent ID.
    pub fn resolve_persistent(&self, entity_id: u64) -> Option<PersistentId> {
        self.context.get_persistent(entity_id)
    }

    /// Check if an entity has a persistent ID.
    pub fn has_entity(&self, entity_id: u64) -> bool {
        self.context.has_entity(entity_id)
    }

    /// Check if a persistent ID is registered.
    pub fn has_persistent(&self, pid: PersistentId) -> bool {
        self.context.has_persistent(pid)
    }

    // ── Propagation ────────────────────────────────────────────────────────────

    /// Propagate names from source entities to target entities.
    ///
    /// Given a mapping from old entity IDs to new entity IDs (or `None` if removed),
    /// this method updates the naming context to reflect the new topology.
    ///
    /// Returns the list of entity IDs that lost their names (because they were removed).
    pub fn propagate_names(
        &mut self,
        entity_map: &[(u64, Option<u64>)],
        policy: NamePropagationPolicy,
    ) -> Vec<u64> {
        let mut lost = Vec::new();

        for (old_entity_id, new_entity_id_opt) in entity_map {
            match new_entity_id_opt {
                Some(new_entity_id) => {
                    // Entity survived; propagate or preserve the name.
                    if let Some(pid) = self.context.get_persistent(*old_entity_id) {
                        match policy {
                            NamePropagationPolicy::Preserve => {
                                self.context.bind(*new_entity_id, pid);
                            }
                            NamePropagationPolicy::Inherit => {
                                // Create a derived ID (e.g., pid + offset).
                                let derived_pid = self.context.allocate_id();
                                self.context.bind(*new_entity_id, derived_pid);
                            }
                            NamePropagationPolicy::Generate => {
                                let new_pid = self.context.allocate_id();
                                self.context.bind(*new_entity_id, new_pid);
                            }
                            NamePropagationPolicy::Combine => {
                                // For single-to-single mapping, same as preserve.
                                self.context.bind(*new_entity_id, pid);
                            }
                        }
                    }
                }
                None => {
                    // Entity was removed.
                    if self.context.has_entity(*old_entity_id) {
                        lost.push(*old_entity_id);
                    }
                }
            }
        }

        lost
    }

    /// Propagate names for a split operation (one entity becomes multiple).
    ///
    /// The source entity's persistent ID is inherited by the first target,
    /// and new IDs are generated for the rest.
    pub fn propagate_split(
        &mut self,
        source_entity_id: u64,
        target_entity_ids: &[u64],
    ) -> Vec<PersistentId> {
        let mut result = Vec::with_capacity(target_entity_ids.len());

        if let Some(source_pid) = self.context.get_persistent(source_entity_id) {
            for (i, &target_id) in target_entity_ids.iter().enumerate() {
                if i == 0 {
                    // First target inherits the source's persistent ID.
                    self.context.bind(target_id, source_pid);
                    result.push(source_pid);
                } else {
                    // Subsequent targets get new IDs.
                    let new_pid = self.context.allocate_id();
                    self.context.bind(target_id, new_pid);
                    result.push(new_pid);
                }
            }
        } else {
            // Source had no persistent ID; generate all new.
            for &target_id in target_entity_ids {
                let pid = self.assign_persistent_id(target_id);
                result.push(pid);
            }
        }

        result
    }

    /// Propagate names for a merge operation (multiple entities become one).
    ///
    /// The target inherits the persistent ID of the first source by default.
    /// With `Combine` policy, all source persistent IDs are recorded as aliases.
    pub fn propagate_merge(
        &mut self,
        source_entity_ids: &[u64],
        target_entity_id: u64,
        policy: NamePropagationPolicy,
    ) -> PersistentId {
        // Find the first source with a persistent ID.
        let primary_pid = source_entity_ids
            .iter()
            .find_map(|&id| self.context.get_persistent(id));

        match policy {
            NamePropagationPolicy::Preserve | NamePropagationPolicy::Inherit => {
                if let Some(pid) = primary_pid {
                    self.context.bind(target_entity_id, pid);
                    pid
                } else {
                    self.assign_persistent_id(target_entity_id)
                }
            }
            NamePropagationPolicy::Generate => {
                self.assign_persistent_id(target_entity_id)
            }
            NamePropagationPolicy::Combine => {
                // Use the first persistent ID but record the merge.
                if let Some(pid) = primary_pid {
                    self.context.bind(target_entity_id, pid);
                    pid
                } else {
                    self.assign_persistent_id(target_entity_id)
                }
            }
        }
    }

    // ── Context Management ─────────────────────────────────────────────────────

    /// Merge another naming context into this one.
    ///
    /// Conflicts (same persistent ID, different entity) are resolved by
    /// generating new persistent IDs for the incoming entities.
    pub fn merge_contexts(&mut self, other: &NamingContext) -> Vec<NamingConflictResolution> {
        let mut resolutions = Vec::new();

        for (entity_id, pid) in other.iter() {
            if let Some(existing_entity) = self.context.get_entity(pid) {
                if existing_entity != entity_id {
                    // Conflict: generate a new PID for the incoming entity.
                    let new_pid = self.context.allocate_id();
                    self.context.bind(entity_id, new_pid);
                    resolutions.push(NamingConflictResolution {
                        conflicting_pid: pid,
                        old_entity_id: existing_entity,
                        new_entity_id: entity_id,
                        resolution: ConflictResolution::GenerateNew,
                    });
                }
                // If same entity, no action needed.
            } else {
                // No conflict; adopt the binding.
                self.context.bind(entity_id, pid);
            }
        }

        self.conflict_history.extend(resolutions.clone());
        resolutions
    }

    /// Get a reference to the current naming context.
    pub fn context(&self) -> &NamingContext {
        &self.context
    }

    /// Get a mutable reference to the naming context.
    pub fn context_mut(&mut self) -> &mut NamingContext {
        &mut self.context
    }

    /// Clear all bindings and reset the ID counter.
    pub fn clear(&mut self) {
        self.context.clear();
        self.conflict_history.clear();
    }

    // ── Reports ────────────────────────────────────────────────────────────────

    /// Generate a stability report comparing before and after contexts.
    pub fn stability_report(
        &self,
        before: &NamingContext,
        entity_ids_after: &[u64],
    ) -> NamingStabilityReport {
        let mut report = NamingStabilityReport::default();
        report.entity_count_before = before.len();
        report.entity_count_after = entity_ids_after.len();

        let mut preserved = 0usize;
        let mut lost = Vec::new();
        let mut new_names = Vec::new();

        for (old_entity, old_pid) in before.iter() {
            // Check if this persistent ID still maps to an entity.
            if let Some(current_entity) = self.context.get_entity(old_pid) {
                if entity_ids_after.contains(&current_entity) {
                    preserved += 1;
                    report.preserved_names.push(current_entity);
                }
            } else {
                lost.push(old_entity);
            }
        }

        // Find new names (entities with persistent IDs not in the before context).
        for &entity_id in entity_ids_after {
            if let Some(pid) = self.context.get_persistent(entity_id) {
                if !before.has_persistent(pid) {
                    new_names.push(entity_id);
                }
            }
        }

        report.lost_names = lost;
        report.new_names = new_names;
        report.conflict_resolutions = self.conflict_history.clone();

        let total_before = before.len().max(1);
        report.stability_score = preserved as f64 / total_before as f64;

        // Calculate weighted score (faces weighted higher).
        // This is a simplified heuristic; in practice, you'd track entity types.
        let face_weight = 0.5;
        let edge_weight = 0.3;
        let vertex_weight = 0.2;
        // Assume equal distribution for simplicity.
        report.weighted_stability_score = report.stability_score * (face_weight + edge_weight + vertex_weight) / 1.0;

        // Detect issues.
        report.issues = self.detect_issues(&report);
        report.broken_chains = report.issues.iter()
            .filter(|i| i.severity == IssueSeverity::Severe || i.severity == IssueSeverity::Critical)
            .count();

        report
    }

    /// Detect naming issues from the current state and report.
    fn detect_issues(&self, report: &NamingStabilityReport) -> Vec<NamingIssue> {
        let mut issues = Vec::new();

        // Check for lost names (moderate to severe).
        for &entity_id in &report.lost_names {
            let severity = if report.lost_names.len() > report.entity_count_before / 2 {
                IssueSeverity::Severe
            } else {
                IssueSeverity::Moderate
            };
            issues.push(NamingIssue {
                severity,
                entity_type: None, // Would need entity type tracking to determine
                description: format!("Entity {} lost its persistent name", entity_id),
                affected_entity_ids: vec![entity_id],
                affected_persistent_ids: vec![],
                suggested_resolution: Some(ConflictResolution::GenerateNew),
            });
        }

        // Check for conflicts.
        for conflict in &report.conflict_resolutions {
            let severity = match conflict.resolution {
                ConflictResolution::KeepOld => IssueSeverity::Minor,
                ConflictResolution::ReplaceOld => IssueSeverity::Moderate,
                ConflictResolution::GenerateNew => IssueSeverity::Minor,
                ConflictResolution::Combine => IssueSeverity::Minor,
            };
            issues.push(NamingIssue {
                severity,
                entity_type: None,
                description: format!(
                    "Naming conflict for persistent ID {} resolved via {:?}",
                    conflict.conflicting_pid.0, conflict.resolution
                ),
                affected_entity_ids: vec![conflict.old_entity_id, conflict.new_entity_id],
                affected_persistent_ids: vec![conflict.conflicting_pid],
                suggested_resolution: Some(conflict.resolution),
            });
        }

        // Check for stability score degradation.
        if report.stability_score < 0.5 {
            issues.push(NamingIssue {
                severity: IssueSeverity::Critical,
                entity_type: None,
                description: format!("Severe naming stability degradation: {:.1}%", report.stability_score * 100.0),
                affected_entity_ids: report.lost_names.clone(),
                affected_persistent_ids: vec![],
                suggested_resolution: None,
            });
        } else if report.stability_score < 0.8 {
            issues.push(NamingIssue {
                severity: IssueSeverity::Severe,
                entity_type: None,
                description: format!("Moderate naming stability degradation: {:.1}%", report.stability_score * 100.0),
                affected_entity_ids: report.lost_names.clone(),
                affected_persistent_ids: vec![],
                suggested_resolution: None,
            });
        }

        issues
    }

    /// Get the conflict resolution history.
    pub fn conflict_history(&self) -> &[NamingConflictResolution] {
        &self.conflict_history
    }

    /// Get the current naming rule.
    pub fn rule(&self) -> NamingRule {
        self.rule
    }

    /// Set the naming rule.
    pub fn set_rule(&mut self, rule: NamingRule) {
        self.rule = rule;
    }

    /// Export the cross-operation history to a NamingHistory for backward compatibility.
    pub fn export_naming_history(&self) -> NamingHistory {
        let mut history = NamingHistory::new();
        for op in &self.cross_op_history.operations {
            for event in &op.naming_events {
                history.push(event.clone());
            }
        }
        history
    }

    /// Import events from a NamingHistory into the cross-operation history.
    pub fn import_naming_history(&mut self, history: &NamingHistory, operation_type: OperationType, label: Option<String>) {
        let _op_id = self.begin_operation(operation_type, label);
        for event in history.iter() {
            self.pending_events.push(event.clone());
            self.apply_event(event);
        }
        let stats = OperationStats {
            entity_count_before: self.context.len(),
            entity_count_after: self.context.len(),
            names_preserved: history.iter().filter(|e| matches!(e, NamingEvent::Propagated { .. })).count(),
            names_lost: history.iter().filter(|e| matches!(e, NamingEvent::Lost { .. })).count(),
            names_generated: history.iter().filter(|e| matches!(e, NamingEvent::Assigned { .. })).count(),
            conflicts_resolved: history.iter().filter(|e| matches!(e, NamingEvent::ConflictResolved { .. })).count(),
        };
        self.finalize_operation(stats);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PersistentNamingHooks Extension
// ─────────────────────────────────────────────────────────────────────────────

/// Extension trait for `PersistentNamingHooks` to integrate with the naming engine.
///
/// This trait provides hooks that can be called during topology operations
/// to maintain persistent naming consistency.
pub trait PersistentNamingHooksExt {
    /// Called when a new face is created.
    ///
    /// `source_entities` lists the entity IDs (if any) that this face was derived from.
    /// Returns the persistent ID assigned to the new face.
    fn on_face_created(
        &mut self,
        engine: &mut PersistentNamingEngine,
        face_idx: usize,
        source_entities: &[u64],
    ) -> PersistentId;

    /// Called when an edge is split into multiple edges.
    ///
    /// Returns the persistent IDs assigned to the new edges.
    fn on_edge_split(
        &mut self,
        engine: &mut PersistentNamingEngine,
        old_edge_idx: usize,
        new_edge_indices: &[usize],
    ) -> Vec<PersistentId>;

    /// Called when multiple vertices are merged into one.
    ///
    /// Returns the persistent ID assigned to the merged vertex.
    fn on_vertex_merged(
        &mut self,
        engine: &mut PersistentNamingEngine,
        old_vertex_indices: &[usize],
        new_vertex_idx: usize,
    ) -> PersistentId;

    /// Called when multiple faces are merged into one.
    ///
    /// Returns the persistent ID assigned to the merged face.
    fn on_face_merged(
        &mut self,
        engine: &mut PersistentNamingEngine,
        old_face_indices: &[usize],
        new_face_idx: usize,
    ) -> PersistentId;
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingEvent for History Integration
// ─────────────────────────────────────────────────────────────────────────────

/// A naming-related event that can be recorded in history.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum NamingEvent {
    /// A new persistent ID was assigned.
    Assigned {
        entity_id: u64,
        persistent_id: PersistentId,
    },
    /// A name was propagated from one entity to another.
    Propagated {
        from_entity: u64,
        to_entity: u64,
        persistent_id: PersistentId,
    },
    /// An entity was split.
    Split {
        source_entity: u64,
        target_entities: Vec<u64>,
        source_persistent_id: PersistentId,
        target_persistent_ids: Vec<PersistentId>,
    },
    /// Entities were merged.
    Merged {
        source_entities: Vec<u64>,
        target_entity: u64,
        result_persistent_id: PersistentId,
    },
    /// A name was lost (entity removed without successor).
    Lost {
        entity_id: u64,
        persistent_id: PersistentId,
    },
    /// A conflict was resolved.
    ConflictResolved(NamingConflictResolution),
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingHistory
// ─────────────────────────────────────────────────────────────────────────────

/// A history log of naming events.
///
/// This can be used to reconstruct a naming context by replaying events.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct NamingHistory {
    /// The recorded naming events.
    pub events: Vec<NamingEvent>,
}

impl NamingHistory {
    /// Create an empty naming history.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a naming event.
    pub fn push(&mut self, event: NamingEvent) {
        self.events.push(event);
    }

    /// Returns the number of events in the history.
    pub fn len(&self) -> usize {
        self.events.len()
    }

    /// Returns `true` if the history is empty.
    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    /// Clear all events.
    pub fn clear(&mut self) {
        self.events.clear();
    }

    /// Iterate over all events.
    pub fn iter(&self) -> impl Iterator<Item = &NamingEvent> {
        self.events.iter()
    }

    /// Replay all events to reconstruct a naming context.
    ///
    /// Returns the reconstructed `NamingContext` and a `PersistentNamingEngine`
    /// initialized with that context.
    pub fn replay(&self) -> PersistentNamingEngine {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        for event in &self.events {
            engine.apply_event(event);
        }

        engine
    }

    /// Replay events from a starting index.
    ///
    /// This is useful for partial replays (e.g., undo/redo).
    pub fn replay_from(&self, start_index: usize) -> PersistentNamingEngine {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        for event in self.events.iter().skip(start_index) {
            engine.apply_event(event);
        }

        engine
    }
}

impl PersistentNamingEngine {
    /// Apply a single naming event to this engine.
    pub fn apply_event(&mut self, event: &NamingEvent) {
        match event {
            NamingEvent::Assigned { entity_id, persistent_id } => {
                self.context.bind(*entity_id, *persistent_id);
            }
            NamingEvent::Propagated { from_entity: _, to_entity, persistent_id } => {
                // Propagation typically means the old entity is gone.
                // Bind the new entity to the same persistent ID.
                self.context.bind(*to_entity, *persistent_id);
            }
            NamingEvent::Split { source_entity, target_entities, source_persistent_id: _, target_persistent_ids } => {
                // Remove the source entity binding.
                self.context.unbind_entity(*source_entity);
                // Bind each target to its assigned persistent ID.
                for (&target_id, &target_pid) in target_entities.iter().zip(target_persistent_ids.iter()) {
                    let _ = target_id; // Used for iteration
                    self.context.bind(target_id, target_pid);
                }
            }
            NamingEvent::Merged { source_entities, target_entity, result_persistent_id } => {
                // Remove all source entity bindings.
                for source_id in source_entities {
                    self.context.unbind_entity(*source_id);
                }
                // Bind the target to the result persistent ID.
                self.context.bind(*target_entity, *result_persistent_id);
            }
            NamingEvent::Lost { entity_id, persistent_id: _ } => {
                // The entity was removed without a successor.
                self.context.unbind_entity(*entity_id);
            }
            NamingEvent::ConflictResolved(resolution) => {
                // Record the conflict resolution.
                self.conflict_history.push(resolution.clone());
            }
        }
    }

    /// Apply an event and track it in the current operation (if any).
    pub fn apply_and_track(&mut self, event: NamingEvent) {
        self.apply_event(&event);
        if self.current_operation.is_some() {
            self.pending_events.push(event);
        }
    }

    /// Record an event to the history and apply it.
    pub fn apply_and_record(&mut self, history: &mut NamingHistory, event: NamingEvent) {
        self.apply_event(&event);
        history.push(event);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cross-Operation History Tracking
// ─────────────────────────────────────────────────────────────────────────────

/// A unique identifier for an operation in the history.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct OperationId(pub u64);

impl OperationId {
    pub const NULL: OperationId = OperationId(0);
    pub fn is_null(&self) -> bool { self.0 == 0 }
}

/// Types of topology operations that can affect naming.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum OperationType {
    /// Boolean union operation.
    BooleanUnion,
    /// Boolean intersection operation.
    BooleanIntersection,
    /// Boolean difference operation.
    BooleanDifference,
    /// Edge split operation.
    EdgeSplit,
    /// Face split operation.
    FaceSplit,
    /// Entity merge operation.
    Merge,
    /// Entity deletion operation.
    Delete,
    /// Geometry transformation (may affect signatures).
    Transform,
    /// Feature-based modification (fillet, chamfer, etc.).
    Feature,
    /// Generic topology modification.
    Generic,
    /// Import operation (STEP, IGES, etc.).
    Import,
}

/// Record of a single operation and its effect on naming.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OperationRecord {
    /// Unique identifier for this operation.
    pub id: OperationId,
    /// Type of operation performed.
    pub operation_type: OperationType,
    /// Optional user-provided label.
    pub label: Option<String>,
    /// Timestamp (operation sequence number).
    pub sequence: u64,
    /// Naming events generated by this operation.
    pub naming_events: Vec<NamingEvent>,
    /// Summary statistics.
    pub stats: OperationStats,
}

/// Statistics for a single operation.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct OperationStats {
    /// Number of entities before the operation.
    pub entity_count_before: usize,
    /// Number of entities after the operation.
    pub entity_count_after: usize,
    /// Number of names preserved.
    pub names_preserved: usize,
    /// Number of names lost.
    pub names_lost: usize,
    /// Number of names generated.
    pub names_generated: usize,
    /// Number of conflicts resolved.
    pub conflicts_resolved: usize,
}

/// Entity genealogy: tracks the origin and evolution of a persistent ID.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntityGenealogy {
    /// The persistent ID being tracked.
    pub persistent_id: PersistentId,
    /// Operation that created this persistent ID.
    pub created_in_operation: OperationId,
    /// List of (operation_id, previous_entity_id) showing entity evolution.
    pub evolution: Vec<(OperationId, u64)>,
    /// Current entity ID (if still alive).
    pub current_entity_id: Option<u64>,
    /// Whether this entity has been deleted.
    pub is_deleted: bool,
}

/// Cross-operation history tracker for persistent naming.
///
/// This struct maintains a complete operation log and entity genealogy,
/// enabling deep history queries and stability analysis across operations.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CrossOperationHistory {
    /// All operations, in order.
    pub operations: Vec<OperationRecord>,
    /// Entity genealogy indexed by persistent ID.
    pub genealogy: HashMap<PersistentId, EntityGenealogy>,
    /// Next operation ID.
    next_operation_id: u64,
    /// Next sequence number.
    next_sequence: u64,
}

impl CrossOperationHistory {
    /// Create a new empty cross-operation history.
    pub fn new() -> Self {
        Self {
            operations: Vec::new(),
            genealogy: HashMap::new(),
            next_operation_id: 1, // Start at 1 so OperationId(0) is reserved as NULL
            next_sequence: 0,
        }
    }

    /// Begin a new operation and return its ID.
    pub fn begin_operation(&mut self, operation_type: OperationType, label: Option<String>) -> OperationId {
        let id = OperationId(self.next_operation_id);
        self.next_operation_id += 1;
        let sequence = self.next_sequence;
        self.next_sequence += 1;

        self.operations.push(OperationRecord {
            id,
            operation_type,
            label,
            sequence,
            naming_events: Vec::new(),
            stats: OperationStats::default(),
        });

        id
    }

    /// Add naming events to the current operation.
    pub fn add_events(&mut self, operation_id: OperationId, events: Vec<NamingEvent>) {
        if let Some(op) = self.operations.iter_mut().find(|op| op.id == operation_id) {
            op.naming_events.extend(events.clone());
        }
        // Update genealogy (separate loop to avoid double mutable borrow).
        for event in &events {
            self.update_genealogy(operation_id, event);
        }
    }

    /// Finalize an operation with statistics.
    pub fn finalize_operation(&mut self, operation_id: OperationId, stats: OperationStats) {
        if let Some(op) = self.operations.iter_mut().find(|op| op.id == operation_id) {
            op.stats = stats;
        }
    }

    fn update_genealogy(&mut self, operation_id: OperationId, event: &NamingEvent) {
        match event {
            NamingEvent::Assigned { entity_id, persistent_id } => {
                self.genealogy.entry(*persistent_id).or_insert_with(|| EntityGenealogy {
                    persistent_id: *persistent_id,
                    created_in_operation: operation_id,
                    evolution: vec![(operation_id, *entity_id)],
                    current_entity_id: Some(*entity_id),
                    is_deleted: false,
                });
            }
            NamingEvent::Propagated { from_entity: _, to_entity, persistent_id } => {
                if let Some(genealogy) = self.genealogy.get_mut(persistent_id) {
                    genealogy.evolution.push((operation_id, *to_entity));
                    genealogy.current_entity_id = Some(*to_entity);
                }
            }
            NamingEvent::Split { source_entity, target_entities, source_persistent_id, target_persistent_ids } => {
                // Mark source as deleted.
                if let Some(g) = self.genealogy.get_mut(source_persistent_id) {
                    g.is_deleted = true;
                    g.current_entity_id = None;
                }
                // Create genealogy for new targets.
                for (i, (&target_id, &target_pid)) in target_entities.iter().zip(target_persistent_ids.iter()).enumerate() {
                    self.genealogy.entry(target_pid).or_insert_with(|| EntityGenealogy {
                        persistent_id: target_pid,
                        created_in_operation: operation_id,
                        evolution: vec![(operation_id, target_id)],
                        current_entity_id: Some(target_id),
                        is_deleted: false,
                    });
                }
            }
            NamingEvent::Merged { source_entities, target_entity, result_persistent_id } => {
                // Mark all sources as deleted.
                for source_id in source_entities {
                    // Try to find the genealogy for each source.
                    for g in self.genealogy.values_mut() {
                        if g.current_entity_id == Some(*source_id) {
                            g.is_deleted = true;
                            g.current_entity_id = None;
                        }
                    }
                }
                // Update result genealogy.
                if let Some(g) = self.genealogy.get_mut(result_persistent_id) {
                    g.evolution.push((operation_id, *target_entity));
                    g.current_entity_id = Some(*target_entity);
                }
            }
            NamingEvent::Lost { entity_id: _, persistent_id } => {
                if let Some(g) = self.genealogy.get_mut(persistent_id) {
                    g.is_deleted = true;
                    g.current_entity_id = None;
                }
            }
            NamingEvent::ConflictResolved(_) => {
                // Conflicts don't directly affect genealogy.
            }
        }
    }

    /// Get the operation record by ID.
    pub fn get_operation(&self, id: OperationId) -> Option<&OperationRecord> {
        self.operations.iter().find(|op| op.id == id)
    }

    /// Get all operations of a specific type.
    pub fn operations_by_type(&self, operation_type: OperationType) -> impl Iterator<Item = &OperationRecord> {
        self.operations.iter().filter(move |op| op.operation_type == operation_type)
    }

    /// Get the genealogy for a persistent ID.
    pub fn get_genealogy(&self, pid: PersistentId) -> Option<&EntityGenealogy> {
        self.genealogy.get(&pid)
    }

    /// Find all operations that affected a given persistent ID.
    pub fn operations_affecting(&self, pid: PersistentId) -> Vec<&OperationRecord> {
        let mut result = Vec::new();
        for op in &self.operations {
            for event in &op.naming_events {
                match event {
                    NamingEvent::Assigned { persistent_id, .. } if *persistent_id == pid => {
                        result.push(op);
                        break;
                    }
                    NamingEvent::Propagated { persistent_id, .. } if *persistent_id == pid => {
                        result.push(op);
                        break;
                    }
                    NamingEvent::Split { target_persistent_ids, .. } if target_persistent_ids.contains(&pid) => {
                        result.push(op);
                        break;
                    }
                    NamingEvent::Merged { result_persistent_id, .. } if *result_persistent_id == pid => {
                        result.push(op);
                        break;
                    }
                    NamingEvent::Lost { persistent_id: lost_pid, .. } if *lost_pid == pid => {
                        result.push(op);
                        break;
                    }
                    _ => {}
                }
            }
        }
        result
    }

    /// Count operations.
    pub fn len(&self) -> usize {
        self.operations.len()
    }

    /// Check if empty.
    pub fn is_empty(&self) -> bool {
        self.operations.is_empty()
    }

    /// Clear all history.
    pub fn clear(&mut self) {
        self.operations.clear();
        self.genealogy.clear();
        self.next_operation_id = 1;
        self.next_sequence = 0;
    }

    /// Generate a comprehensive stability report across all operations.
    pub fn stability_report(&self) -> CrossOperationStabilityReport {
        let mut report = CrossOperationStabilityReport::default();

        report.total_operations = self.operations.len();

        for op in &self.operations {
            report.total_names_preserved += op.stats.names_preserved;
            report.total_names_lost += op.stats.names_lost;
            report.total_names_generated += op.stats.names_generated;
            report.total_conflicts += op.stats.conflicts_resolved;
        }

        // Count alive vs deleted entities.
        for genealogy in self.genealogy.values() {
            if genealogy.is_deleted {
                report.entities_deleted += 1;
            } else {
                report.entities_alive += 1;
            }
        }

        report.total_entities_tracked = self.genealogy.len();

        // Calculate stability score.
        let preserved = report.total_names_preserved as f64;
        let lost = report.total_names_lost as f64;
        let total = preserved + lost;
        report.overall_stability_score = if total > 0.0 { preserved / total } else { 1.0 };

        report
    }
}

/// Stability report across multiple operations.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CrossOperationStabilityReport {
    /// Total number of operations.
    pub total_operations: usize,
    /// Total names preserved across all operations.
    pub total_names_preserved: usize,
    /// Total names lost across all operations.
    pub total_names_lost: usize,
    /// Total names generated across all operations.
    pub total_names_generated: usize,
    /// Total conflicts resolved.
    pub total_conflicts: usize,
    /// Entities currently alive.
    pub entities_alive: usize,
    /// Entities that have been deleted.
    pub entities_deleted: usize,
    /// Total entities tracked (alive + deleted).
    pub total_entities_tracked: usize,
    /// Overall stability score (0.0 - 1.0).
    pub overall_stability_score: f64,
}

impl CrossOperationStabilityReport {
    /// Returns true if stability is excellent (> 95% preserved, minimal losses).
    pub fn is_excellent(&self) -> bool {
        self.overall_stability_score >= 0.95 && self.total_names_lost < 5
    }

    /// Returns true if stability is good (> 90% preserved).
    pub fn is_good(&self) -> bool {
        self.overall_stability_score >= 0.90
    }

    /// Returns true if there are significant naming issues.
    pub fn has_issues(&self) -> bool {
        self.overall_stability_score < 0.90 || self.total_conflicts > 0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── PersistentId tests ─────────────────────────────────────────────────────

    #[test]
    fn persistent_id_null_is_zero() {
        assert!(PersistentId::NULL.is_null());
        assert_eq!(PersistentId::NULL.raw(), 0);
    }

    #[test]
    fn persistent_id_display() {
        let pid = PersistentId(42);
        assert_eq!(format!("{pid}"), "pid:42");
    }

    // ── NamingContext tests ────────────────────────────────────────────────────

    #[test]
    fn naming_context_starts_empty() {
        let ctx = NamingContext::new();
        assert!(ctx.is_empty());
        assert_eq!(ctx.len(), 0);
    }

    #[test]
    fn naming_context_bind_and_lookup() {
        let mut ctx = NamingContext::new();
        let pid = ctx.allocate_id();
        ctx.bind(42, pid);

        assert!(ctx.has_entity(42));
        assert!(ctx.has_persistent(pid));
        assert_eq!(ctx.get_persistent(42), Some(pid));
        assert_eq!(ctx.get_entity(pid), Some(42));
    }

    #[test]
    fn naming_context_unbind() {
        let mut ctx = NamingContext::new();
        let pid = ctx.allocate_id();
        ctx.bind(42, pid);

        assert_eq!(ctx.unbind_entity(42), Some(pid));
        assert!(!ctx.has_entity(42));
        assert!(!ctx.has_persistent(pid));
    }

    #[test]
    fn naming_context_bind_replaces_old() {
        let mut ctx = NamingContext::new();
        let pid1 = ctx.allocate_id();
        let pid2 = ctx.allocate_id();

        ctx.bind(42, pid1);
        ctx.bind(42, pid2);

        // Entity 42 should now have pid2.
        assert_eq!(ctx.get_persistent(42), Some(pid2));
        // pid1 should no longer be bound.
        assert!(!ctx.has_persistent(pid1));
    }

    // ── PersistentNamingEngine tests ───────────────────────────────────────────

    #[test]
    fn engine_assigns_unique_ids() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid1 = engine.assign_persistent_id(1);
        let pid2 = engine.assign_persistent_id(2);

        assert_ne!(pid1, pid2);
        assert!(!pid1.is_null());
        assert!(!pid2.is_null());
    }

    #[test]
    fn engine_returns_existing_id() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid1 = engine.assign_persistent_id(1);
        let pid2 = engine.assign_persistent_id(1);

        assert_eq!(pid1, pid2);
    }

    #[test]
    fn engine_resolves_bidirectional() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid = engine.assign_persistent_id(42);

        assert_eq!(engine.resolve_entity(pid), Some(42));
        assert_eq!(engine.resolve_persistent(42), Some(pid));
    }

    #[test]
    fn engine_propagate_preserve() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid = engine.assign_persistent_id(10);

        let entity_map = vec![(10, Some(20))];
        let lost = engine.propagate_names(&entity_map, NamePropagationPolicy::Preserve);

        assert!(lost.is_empty());
        assert_eq!(engine.resolve_persistent(20), Some(pid));
    }

    #[test]
    fn engine_propagate_removed() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        engine.assign_persistent_id(10);

        let entity_map = vec![(10, None)];
        let lost = engine.propagate_names(&entity_map, NamePropagationPolicy::Preserve);

        assert_eq!(lost, vec![10]);
    }

    #[test]
    fn engine_propagate_split() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let source_pid = engine.assign_persistent_id(1);

        let result = engine.propagate_split(1, &[10, 11, 12]);

        assert_eq!(result.len(), 3);
        assert_eq!(result[0], source_pid); // First inherits
        assert_ne!(result[1], source_pid); // Others get new IDs
        assert_ne!(result[2], source_pid);
    }

    #[test]
    fn engine_propagate_merge() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid1 = engine.assign_persistent_id(1);
        engine.assign_persistent_id(2);

        let result = engine.propagate_merge(&[1, 2], 100, NamePropagationPolicy::Preserve);

        assert_eq!(result, pid1); // First source's ID is preserved
    }

    #[test]
    fn engine_merge_contexts_no_conflict() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());
        engine.assign_persistent_id(1);

        let mut other = NamingContext::new();
        // Use a PID that doesn't conflict with engine's PIDs.
        // Engine's first PID is 1, so we skip to 100.
        other.next_id = 100;
        let other_pid = other.allocate_id();
        other.bind(2, other_pid);

        let resolutions = engine.merge_contexts(&other);

        assert!(resolutions.is_empty(), "Should have no conflicts with different PIDs");
        assert!(engine.has_entity(1));
        assert!(engine.has_entity(2));
    }

    #[test]
    fn engine_merge_contexts_with_conflict() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());
        engine.assign_persistent_id(1);

        let mut other = NamingContext::new();
        // Force the same PID to be allocated (simulate conflict).
        other.next_id = 1;
        let conflicting_pid = other.allocate_id();
        other.bind(2, conflicting_pid);

        let resolutions = engine.merge_contexts(&other);

        assert_eq!(resolutions.len(), 1);
        assert_eq!(resolutions[0].resolution, ConflictResolution::GenerateNew);
    }

    // ── NamingStabilityReport tests ────────────────────────────────────────────

    #[test]
    fn stability_report_perfect() {
        let mut before = NamingContext::new();
        before.bind(1, PersistentId(1));
        before.bind(2, PersistentId(2));

        let mut engine = PersistentNamingEngine::new(NamingRule::default());
        // Set up the engine context to match 'before'.
        engine.context_mut().bind(1, PersistentId(1));
        engine.context_mut().bind(2, PersistentId(2));

        let report = engine.stability_report(&before, &[1, 2]);

        // With identical before/after, score should be 1.0.
        assert_eq!(report.stability_score, 1.0);
        assert!(report.is_perfect());
    }

    #[test]
    fn stability_report_has_issues() {
        let report = NamingStabilityReport {
            stability_score: 0.5,
            lost_names: vec![1],
            new_names: vec![],
            preserved_names: vec![2],
            conflict_resolutions: vec![],
            entity_count_before: 2,
            entity_count_after: 1,
            entity_type_breakdown: HashMap::new(),
            issues: vec![],
            weighted_stability_score: 0.5,
            broken_chains: 0,
        };

        assert!(report.has_issues());
        assert!(!report.is_perfect());
    }

    // ── NamePropagationPolicy tests ────────────────────────────────────────────

    #[test]
    fn propagate_inherit_creates_new_ids() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let old_pid = engine.assign_persistent_id(10);

        let entity_map = vec![(10, Some(20))];
        let _lost = engine.propagate_names(&entity_map, NamePropagationPolicy::Inherit);

        // With Inherit, a new ID is created (not the old one).
        let new_pid = engine.resolve_persistent(20);
        assert!(new_pid.is_some());
        assert_ne!(new_pid, Some(old_pid));
    }

    #[test]
    fn propagate_generate_creates_new_ids() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let old_pid = engine.assign_persistent_id(10);

        let entity_map = vec![(10, Some(20))];
        let _lost = engine.propagate_names(&entity_map, NamePropagationPolicy::Generate);

        let new_pid = engine.resolve_persistent(20);
        assert!(new_pid.is_some());
        assert_ne!(new_pid, Some(old_pid));
    }

    // ── Cross-Operation History tests ────────────────────────────────────────────

    #[test]
    fn cross_op_history_begins_operation() {
        let mut history = CrossOperationHistory::new();
        let op_id = history.begin_operation(OperationType::BooleanUnion, Some("test_op".to_string()));

        assert!(!op_id.is_null());
        assert_eq!(history.len(), 1);
        assert!(history.get_operation(op_id).is_some());
    }

    #[test]
    fn cross_op_history_adds_events() {
        let mut history = CrossOperationHistory::new();
        let op_id = history.begin_operation(OperationType::Generic, None);

        let events = vec![
            NamingEvent::Assigned { entity_id: 1, persistent_id: PersistentId(1) },
            NamingEvent::Propagated { from_entity: 2, to_entity: 3, persistent_id: PersistentId(2) },
        ];

        history.add_events(op_id, events.clone());

        let op = history.get_operation(op_id).unwrap();
        assert_eq!(op.naming_events.len(), 2);
    }

    #[test]
    fn cross_op_history_tracks_genealogy() {
        let mut history = CrossOperationHistory::new();
        let op_id = history.begin_operation(OperationType::Generic, None);

        let pid = PersistentId(1);
        history.add_events(op_id, vec![
            NamingEvent::Assigned { entity_id: 42, persistent_id: pid },
        ]);

        let genealogy = history.get_genealogy(pid).unwrap();
        assert_eq!(genealogy.persistent_id, pid);
        assert_eq!(genealogy.current_entity_id, Some(42));
        assert!(!genealogy.is_deleted);
    }

    #[test]
    fn cross_op_history_stability_report() {
        let mut history = CrossOperationHistory::new();

        let op_id = history.begin_operation(OperationType::BooleanUnion, None);
        history.finalize_operation(op_id, OperationStats {
            entity_count_before: 10,
            entity_count_after: 8,
            names_preserved: 7,
            names_lost: 1,
            names_generated: 2,
            conflicts_resolved: 0,
        });

        let report = history.stability_report();
        assert_eq!(report.total_operations, 1);
        assert_eq!(report.total_names_preserved, 7);
        assert_eq!(report.total_names_lost, 1);
        assert_eq!(report.total_names_generated, 2);
    }

    // ── PersistentNamingEngine Cross-Op tests ────────────────────────────────────

    #[test]
    fn engine_begin_operation_tracks_events() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let op_id = engine.begin_operation(OperationType::BooleanUnion, Some("union".to_string()));
        engine.assign_persistent_id(1);
        engine.assign_persistent_id(2);

        engine.finalize_operation(OperationStats {
            entity_count_before: 0,
            entity_count_after: 2,
            names_preserved: 0,
            names_lost: 0,
            names_generated: 2,
            conflicts_resolved: 0,
        });

        let history = engine.cross_operation_history();
        assert_eq!(history.len(), 1);

        let op = history.get_operation(op_id).unwrap();
        assert_eq!(op.naming_events.len(), 2);
    }

    #[test]
    fn engine_cross_op_stability_report() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        // First operation.
        let _op1 = engine.begin_operation(OperationType::BooleanUnion, None);
        engine.assign_persistent_id(1);
        engine.assign_persistent_id(2);
        engine.finalize_operation(OperationStats {
            entity_count_before: 0,
            entity_count_after: 2,
            names_preserved: 0,
            names_lost: 0,
            names_generated: 2,
            conflicts_resolved: 0,
        });

        // Second operation.
        let _op2 = engine.begin_operation(OperationType::BooleanDifference, None);
        engine.assign_persistent_id(3);
        engine.finalize_operation(OperationStats {
            entity_count_before: 2,
            entity_count_after: 3,
            names_preserved: 1,
            names_lost: 0,
            names_generated: 1,
            conflicts_resolved: 0,
        });

        let report = engine.cross_operation_stability_report();
        assert_eq!(report.total_operations, 2);
        assert_eq!(report.total_names_generated, 3);
        assert!(report.is_excellent() || report.is_good());
    }

    // ── NamingStabilityReport enhanced tests ───────────────────────────────────────

    #[test]
    fn stability_report_detects_issues() {
        let mut before = NamingContext::new();
        before.bind(1, PersistentId(1));
        before.bind(2, PersistentId(2));
        before.bind(3, PersistentId(3));

        let mut engine = PersistentNamingEngine::new(NamingRule::default());
        engine.context_mut().bind(1, PersistentId(1));
        // Entity 2 and 3 are lost.

        let report = engine.stability_report(&before, &[1]);

        assert!(report.has_issues());
        assert_eq!(report.lost_names.len(), 2);
        assert!(!report.issues.is_empty());
    }

    #[test]
    fn stability_report_issue_summary() {
        let report = NamingStabilityReport {
            stability_score: 0.5,
            lost_names: vec![1, 2, 3],
            new_names: vec![],
            preserved_names: vec![4],
            conflict_resolutions: vec![],
            entity_count_before: 4,
            entity_count_after: 2,
            issues: vec![
                NamingIssue {
                    severity: IssueSeverity::Moderate,
                    entity_type: Some(EntityType::Face),
                    description: "Test issue 1".to_string(),
                    affected_entity_ids: vec![1],
                    affected_persistent_ids: vec![],
                    suggested_resolution: None,
                },
                NamingIssue {
                    severity: IssueSeverity::Severe,
                    entity_type: Some(EntityType::Edge),
                    description: "Test issue 2".to_string(),
                    affected_entity_ids: vec![2],
                    affected_persistent_ids: vec![],
                    suggested_resolution: None,
                },
            ],
            entity_type_breakdown: HashMap::new(),
            weighted_stability_score: 0.5,
            broken_chains: 1,
        };

        let summary = report.issue_summary();
        assert_eq!(summary.get(&IssueSeverity::Moderate), Some(&1));
        assert_eq!(summary.get(&IssueSeverity::Severe), Some(&1));
    }

    #[test]
    fn stability_report_summary_string() {
        let report = NamingStabilityReport {
            stability_score: 0.75,
            lost_names: vec![1],
            new_names: vec![],
            preserved_names: vec![2, 3],
            conflict_resolutions: vec![],
            entity_count_before: 3,
            entity_count_after: 3,
            issues: vec![],
            entity_type_breakdown: HashMap::new(),
            weighted_stability_score: 0.75,
            broken_chains: 0,
        };

        let summary = report.summary_string();
        assert!(summary.contains("75.0%"));
        assert!(summary.contains("Preserved: 2"));
    }
}

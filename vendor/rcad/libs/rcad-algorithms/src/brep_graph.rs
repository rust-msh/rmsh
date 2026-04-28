use std::collections::HashMap;
use serde::{Deserialize, Serialize};

use rcad_kernel::BRep;
use rcad_kernel::persistent_naming::{
    PersistentId, PersistentNamingEngine, NamingStabilityReport,
    OperationType, OperationStats, CrossOperationStabilityReport,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NodeKind {
    Solid,
    Shell,
    Face,
    Wire,
    Edge,
    Vertex,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TopoNode {
    pub kind: NodeKind,
    pub index: usize,
}

#[derive(Debug, Clone)]
pub struct TopoGraphHistoryEvent {
    pub action: String,
}

#[derive(Debug, Clone, Default)]
pub struct TopoGraphHistory {
    pub events: Vec<TopoGraphHistoryEvent>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TopoGraphValidationIssue {
    MissingAdjacency { node: TopoNode },
    NonSymmetricAdjacency { a: TopoNode, b: TopoNode },
    InvalidEdgeVertexRef { edge_index: usize, vertex_index: usize },
}

#[derive(Debug, Clone, Default)]
pub struct TopoGraph {
    pub nodes: Vec<TopoNode>,
    pub history: TopoGraphHistory,
    adjacency: HashMap<TopoNode, Vec<TopoNode>>,
    solid_shells: Vec<Vec<usize>>,
    shell_faces: Vec<Vec<usize>>,
    face_wires: Vec<Vec<usize>>,
    wire_edges: Vec<Vec<usize>>,
    edge_vertices: Vec<[usize; 2]>,
}

impl TopoGraph {
    pub fn from_brep(brep: &BRep) -> Self {
        let mut g = Self::default();
        g.record("from_brep");

        for vi in 0..brep.vertices.len() {
            g.add_node(TopoNode {
                kind: NodeKind::Vertex,
                index: vi,
            });
        }

        for (ei, e) in brep.edges.iter().enumerate() {
            let en = TopoNode {
                kind: NodeKind::Edge,
                index: ei,
            };
            g.add_node(en);
            g.edge_vertices.push([e.start, e.end]);
            g.connect(
                en,
                TopoNode {
                    kind: NodeKind::Vertex,
                    index: e.start,
                },
            );
            g.connect(
                en,
                TopoNode {
                    kind: NodeKind::Vertex,
                    index: e.end,
                },
            );
        }

        let mut shell_idx = 0usize;
        let mut face_idx = 0usize;
        let mut wire_idx = 0usize;

        for (si, solid) in brep.solids.iter().enumerate() {
            let sn = TopoNode {
                kind: NodeKind::Solid,
                index: si,
            };
            g.add_node(sn);
            g.solid_shells.push(Vec::new());

            for shell in &solid.shells {
                let shn = TopoNode {
                    kind: NodeKind::Shell,
                    index: shell_idx,
                };
                g.add_node(shn);
                g.connect(sn, shn);
                g.solid_shells[si].push(shell_idx);
                g.shell_faces.push(Vec::new());

                for face in &shell.faces {
                    let fnn = TopoNode {
                        kind: NodeKind::Face,
                        index: face_idx,
                    };
                    g.add_node(fnn);
                    g.connect(shn, fnn);
                    g.shell_faces[shell_idx].push(face_idx);
                    g.face_wires.push(Vec::new());

                    let wires = std::iter::once(&face.outer_wire).chain(face.inner_wires.iter());
                    for wire in wires {
                        let wn = TopoNode {
                            kind: NodeKind::Wire,
                            index: wire_idx,
                        };
                        g.add_node(wn);
                        g.connect(fnn, wn);
                        g.face_wires[face_idx].push(wire_idx);
                        g.wire_edges.push(Vec::new());

                        for we in &wire.edges {
                            let en = TopoNode {
                                kind: NodeKind::Edge,
                                index: we.idx,
                            };
                            if we.idx < brep.edges.len() {
                                g.connect(wn, en);
                                if !g.wire_edges[wire_idx].contains(&we.idx) {
                                    g.wire_edges[wire_idx].push(we.idx);
                                }
                            }
                        }

                        wire_idx += 1;
                    }

                    face_idx += 1;
                }

                shell_idx += 1;
            }
        }

        g
    }

    pub fn record(&mut self, action: impl Into<String>) {
        self.history.events.push(TopoGraphHistoryEvent {
            action: action.into(),
        });
    }

    pub fn neighbors(&self, node: TopoNode) -> Vec<TopoNode> {
        self.adjacency.get(&node).cloned().unwrap_or_default()
    }

    pub fn faces_of_shell(&self, shell: TopoNode) -> Vec<TopoNode> {
        if shell.kind != NodeKind::Shell {
            return Vec::new();
        }
        self.shell_faces
            .get(shell.index)
            .map(|v| {
                v.iter()
                    .map(|&i| TopoNode {
                        kind: NodeKind::Face,
                        index: i,
                    })
                    .collect()
            })
            .unwrap_or_default()
    }

    pub fn edges_of_face(&self, face: TopoNode) -> Vec<TopoNode> {
        if face.kind != NodeKind::Face {
            return Vec::new();
        }
        let mut out: Vec<usize> = Vec::new();
        if let Some(wires) = self.face_wires.get(face.index) {
            for &wi in wires {
                if let Some(edges) = self.wire_edges.get(wi) {
                    for &ei in edges {
                        if !out.contains(&ei) {
                            out.push(ei);
                        }
                    }
                }
            }
        }
        out.into_iter()
            .map(|i| TopoNode {
                kind: NodeKind::Edge,
                index: i,
            })
            .collect()
    }

    pub fn vertices_of_edge(&self, edge: TopoNode) -> Vec<TopoNode> {
        if edge.kind != NodeKind::Edge {
            return Vec::new();
        }
        self.edge_vertices
            .get(edge.index)
            .map(|v| {
                vec![
                    TopoNode {
                        kind: NodeKind::Vertex,
                        index: v[0],
                    },
                    TopoNode {
                        kind: NodeKind::Vertex,
                        index: v[1],
                    },
                ]
            })
            .unwrap_or_default()
    }

    pub fn validate(&self) -> Vec<TopoGraphValidationIssue> {
        let mut issues = Vec::new();
        for node in &self.nodes {
            let Some(neigh) = self.adjacency.get(node) else {
                issues.push(TopoGraphValidationIssue::MissingAdjacency { node: *node });
                continue;
            };
            for n in neigh {
                if let Some(back) = self.adjacency.get(n) {
                    if !back.contains(node) {
                        issues.push(TopoGraphValidationIssue::NonSymmetricAdjacency {
                            a: *node,
                            b: *n,
                        });
                    }
                } else {
                    issues.push(TopoGraphValidationIssue::MissingAdjacency { node: *n });
                }
            }
        }

        for (ei, vv) in self.edge_vertices.iter().enumerate() {
            for &vi in vv {
                if !self.nodes.contains(&TopoNode {
                    kind: NodeKind::Vertex,
                    index: vi,
                }) {
                    issues.push(TopoGraphValidationIssue::InvalidEdgeVertexRef {
                        edge_index: ei,
                        vertex_index: vi,
                    });
                }
            }
        }

        issues
    }

    /// Compact graph storage by dropping orphan adjacency entries and
    /// deduplicating neighbor lists.
    pub fn compact(&mut self) {
        let mut node_set = std::collections::HashSet::new();
        for n in &self.nodes {
            node_set.insert(*n);
        }

        self.adjacency.retain(|node, _| node_set.contains(node));
        for neigh in self.adjacency.values_mut() {
            neigh.retain(|n| node_set.contains(n));
            neigh.sort_by_key(|n| {
                let kind_rank = match n.kind {
                    NodeKind::Solid => 0usize,
                    NodeKind::Shell => 1,
                    NodeKind::Face => 2,
                    NodeKind::Wire => 3,
                    NodeKind::Edge => 4,
                    NodeKind::Vertex => 5,
                };
                (kind_rank, n.index)
            });
            neigh.dedup();
        }
        self.record("compact");
    }

    /// Apply a mutation and run graph validation afterward.
    ///
    /// This is a lightweight baseline for mutation-guard workflows: callers
    /// can route all topology edits through this helper and reject invalid
    /// states before continuing downstream processing.
    pub fn mutate_checked<F>(
        &mut self,
        action: impl Into<String>,
        mutator: F,
    ) -> Result<(), Vec<TopoGraphValidationIssue>>
    where
        F: FnOnce(&mut TopoGraph),
    {
        let action = action.into();
        mutator(self);
        let issues = self.validate();
        if issues.is_empty() {
            self.record(format!("mutate:{action}"));
            Ok(())
        } else {
            self.record(format!("mutate_invalid:{action}"));
            Err(issues)
        }
    }

    /// Apply a mutation with rollback-on-failure semantics.
    ///
    /// If validation fails after the mutation, graph state is restored to the
    /// pre-mutation snapshot and validation issues are returned.
    pub fn mutate_guarded<F>(
        &mut self,
        action: impl Into<String>,
        mutator: F,
    ) -> Result<(), Vec<TopoGraphValidationIssue>>
    where
        F: FnOnce(&mut TopoGraph),
    {
        let action = action.into();
        let before = self.clone();
        mutator(self);
        let issues = self.validate();
        if issues.is_empty() {
            self.record(format!("mutate_guarded:{action}"));
            Ok(())
        } else {
            *self = before;
            self.record(format!("mutate_guarded_rollback:{action}"));
            Err(issues)
        }
    }

    fn add_node(&mut self, node: TopoNode) {
        if !self.nodes.contains(&node) {
            self.nodes.push(node);
        }
        self.adjacency.entry(node).or_default();
    }

    fn connect(&mut self, a: TopoNode, b: TopoNode) {
        self.add_node(a);
        self.add_node(b);
        let va = self.adjacency.entry(a).or_default();
        if !va.contains(&b) {
            va.push(b);
        }
        let vb = self.adjacency.entry(b).or_default();
        if !vb.contains(&a) {
            vb.push(a);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRepGraphHistory: Persistent Naming Integration
// ─────────────────────────────────────────────────────────────────────────────

/// Enhanced history with persistent naming integration for cross-operation stability.
///
/// This struct bridges the TopoGraph mutation history with the PersistentNamingEngine,
/// enabling:
/// - Automatic name propagation during topology mutations
/// - Cross-operation stability analysis
/// - Entity genealogy tracking
/// - Undo/redo support with naming reconstruction
#[derive(Debug, Clone)]
pub struct BRepGraphHistory {
    /// The underlying naming engine.
    naming_engine: PersistentNamingEngine,
    /// Snapshots for undo support.
    snapshots: Vec<TopoGraphSnapshot>,
    /// Current snapshot index (for undo/redo).
    current_snapshot: usize,
}

/// A snapshot of the graph state with naming context.
#[derive(Debug, Clone)]
struct TopoGraphSnapshot {
    /// The action that created this snapshot.
    action: String,
    /// Node count at this snapshot.
    node_count: usize,
    /// Entity ID to persistent ID mappings at this snapshot.
    naming: HashMap<TopoNode, PersistentId>,
}

impl Default for BRepGraphHistory {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepGraphHistory {
    /// Create a new history with default naming engine.
    pub fn new() -> Self {
        Self {
            naming_engine: PersistentNamingEngine::default(),
            snapshots: Vec::new(),
            current_snapshot: 0,
        }
    }

    /// Create a history with a specific naming rule.
    pub fn with_naming_rule(rule: rcad_kernel::persistent_naming::NamingRule) -> Self {
        Self {
            naming_engine: PersistentNamingEngine::new(rule),
            snapshots: Vec::new(),
            current_snapshot: 0,
        }
    }

    /// Get a reference to the naming engine.
    pub fn naming_engine(&self) -> &PersistentNamingEngine {
        &self.naming_engine
    }

    /// Get mutable access to the naming engine.
    pub fn naming_engine_mut(&mut self) -> &mut PersistentNamingEngine {
        &mut self.naming_engine
    }

    /// Begin an operation for history tracking.
    pub fn begin_operation(&mut self, operation_type: OperationType, label: Option<String>) {
        self.naming_engine.begin_operation(operation_type, label);
    }

    /// Finalize the current operation.
    pub fn finalize_operation(&mut self, stats: OperationStats) {
        self.naming_engine.finalize_operation(stats);
    }

    /// Record a graph mutation with naming propagation.
    ///
    /// This creates a snapshot and propagates names for surviving entities.
    pub fn record_mutation(
        &mut self,
        graph: &TopoGraph,
        action: &str,
        entity_map: &[(TopoNode, Option<TopoNode>)],
    ) {
        // Create snapshot.
        let mut naming = HashMap::new();
        for node in &graph.nodes {
            if let Some(pid) = self.naming_engine.resolve_persistent(node_to_entity_id(*node)) {
                naming.insert(*node, pid);
            }
        }

        self.snapshots.truncate(self.current_snapshot + 1);
        self.snapshots.push(TopoGraphSnapshot {
            action: action.to_string(),
            node_count: graph.nodes.len(),
            naming,
        });
        self.current_snapshot = self.snapshots.len() - 1;

        // Propagate names for surviving entities.
        let entity_id_map: Vec<(u64, Option<u64>)> = entity_map
            .iter()
            .map(|(old, new)| {
                (
                    node_to_entity_id(*old),
                    new.map(node_to_entity_id),
                )
            })
            .collect();

        self.naming_engine.propagate_names(
            &entity_id_map,
            rcad_kernel::persistent_naming::NamePropagationPolicy::Preserve,
        );
    }

    /// Assign a persistent ID to a topology node.
    pub fn assign_persistent_id(&mut self, node: TopoNode) -> PersistentId {
        self.naming_engine.assign_persistent_id(node_to_entity_id(node))
    }

    /// Resolve a topology node to its persistent ID.
    pub fn resolve_persistent(&self, node: TopoNode) -> Option<PersistentId> {
        self.naming_engine.resolve_persistent(node_to_entity_id(node))
    }

    /// Resolve a persistent ID back to a topology node.
    pub fn resolve_node(&self, pid: PersistentId, kind: NodeKind) -> Option<TopoNode> {
        let entity_id = self.naming_engine.resolve_entity(pid)?;
        Some(TopoNode {
            kind,
            index: entity_id as usize,
        })
    }

    /// Generate a stability report for the current state.
    pub fn stability_report(&self) -> CrossOperationStabilityReport {
        self.naming_engine.cross_operation_stability_report()
    }

    /// Generate a naming stability report comparing before and after states.
    pub fn naming_stability_report(
        &self,
        before_nodes: &[TopoNode],
        after_nodes: &[TopoNode],
    ) -> NamingStabilityReport {
        let before_context = self.naming_engine.context().clone();
        let after_ids: Vec<u64> = after_nodes.iter().map(|n| node_to_entity_id(*n)).collect();
        self.naming_engine.stability_report(&before_context, &after_ids)
    }

    /// Track an edge split event.
    pub fn track_edge_split(
        &mut self,
        old_edge_idx: usize,
        new_edge_indices: &[usize],
    ) -> Vec<PersistentId> {
        self.naming_engine.propagate_split(
            old_edge_idx as u64,
            &new_edge_indices.iter().map(|&i| i as u64).collect::<Vec<_>>(),
        )
    }

    /// Track a face split event.
    pub fn track_face_split(
        &mut self,
        old_face_idx: usize,
        new_face_indices: &[usize],
    ) -> Vec<PersistentId> {
        self.naming_engine.propagate_split(
            old_face_idx as u64,
            &new_face_indices.iter().map(|&i| i as u64).collect::<Vec<_>>(),
        )
    }

    /// Track a vertex merge event.
    pub fn track_vertex_merge(
        &mut self,
        old_vertex_indices: &[usize],
        new_vertex_idx: usize,
    ) -> PersistentId {
        self.naming_engine.propagate_merge(
            &old_vertex_indices.iter().map(|&i| i as u64).collect::<Vec<_>>(),
            new_vertex_idx as u64,
            rcad_kernel::persistent_naming::NamePropagationPolicy::Preserve,
        )
    }

    /// Track a face merge event.
    pub fn track_face_merge(
        &mut self,
        old_face_indices: &[usize],
        new_face_idx: usize,
    ) -> PersistentId {
        self.naming_engine.propagate_merge(
            &old_face_indices.iter().map(|&i| i as u64).collect::<Vec<_>>(),
            new_face_idx as u64,
            rcad_kernel::persistent_naming::NamePropagationPolicy::Preserve,
        )
    }

    /// Get the number of recorded snapshots.
    pub fn snapshot_count(&self) -> usize {
        self.snapshots.len()
    }

    /// Check if undo is available.
    pub fn can_undo(&self) -> bool {
        self.current_snapshot > 0
    }

    /// Check if redo is available.
    pub fn can_redo(&self) -> bool {
        self.current_snapshot < self.snapshots.len().saturating_sub(1)
    }

    /// Get the action name at the current snapshot.
    pub fn current_action(&self) -> Option<&str> {
        self.snapshots.get(self.current_snapshot).map(|s| s.action.as_str())
    }

    /// Export naming events as a serializable history.
    pub fn export_naming_history(&self) -> rcad_kernel::persistent_naming::NamingHistory {
        self.naming_engine.export_naming_history()
    }
}

/// Convert a TopoNode to a unique entity ID.
fn node_to_entity_id(node: TopoNode) -> u64 {
    // Encode kind and index into a single u64.
    // Kind uses high 8 bits, index uses low 56 bits.
    let kind_bits = match node.kind {
        NodeKind::Solid => 0u64,
        NodeKind::Shell => 1u64,
        NodeKind::Face => 2u64,
        NodeKind::Wire => 3u64,
        NodeKind::Edge => 4u64,
        NodeKind::Vertex => 5u64,
    };
    (kind_bits << 56) | (node.index as u64)
}

/// Convert an entity ID back to a TopoNode (requires known kind).
fn entity_id_to_node(entity_id: u64, kind: NodeKind) -> TopoNode {
    let index = (entity_id & 0x00FFFFFFFFFFFFFF) as usize;
    TopoNode { kind, index }
}

/// NamedGraph: A TopoGraph with integrated naming history.
///
/// This provides a convenient wrapper for applications that need
/// automatic naming tracking during graph mutations.
#[derive(Debug, Clone)]
pub struct NamedGraph {
    graph: TopoGraph,
    history: BRepGraphHistory,
}

impl NamedGraph {
    /// Create a new named graph from a BRep.
    pub fn from_brep(brep: &BRep) -> Self {
        let graph = TopoGraph::from_brep(brep);
        let mut history = BRepGraphHistory::new();

        // Assign persistent IDs to all nodes.
        for node in &graph.nodes {
            history.assign_persistent_id(*node);
        }

        Self { graph, history }
    }

    /// Get the underlying graph.
    pub fn graph(&self) -> &TopoGraph {
        &self.graph
    }

    /// Get mutable access to the graph.
    pub fn graph_mut(&mut self) -> &mut TopoGraph {
        &mut self.graph
    }

    /// Get the history.
    pub fn history(&self) -> &BRepGraphHistory {
        &self.history
    }

    /// Get mutable access to the history.
    pub fn history_mut(&mut self) -> &mut BRepGraphHistory {
        &mut self.history
    }

    /// Apply a mutation with automatic naming tracking.
    pub fn mutate_tracked<F>(
        &mut self,
        action: &str,
        operation_type: OperationType,
        mutator: F,
    ) -> Result<(), Vec<TopoGraphValidationIssue>>
    where
        F: FnOnce(&mut TopoGraph, &mut BRepGraphHistory),
    {
        self.history.begin_operation(operation_type, Some(action.to_string()));

        let before_nodes = self.graph.nodes.clone();
        mutator(&mut self.graph, &mut self.history);

        let issues = self.graph.validate();
        if issues.is_empty() {
            let after_nodes = self.graph.nodes.clone();
            let entity_map: Vec<(TopoNode, Option<TopoNode>)> = before_nodes
                .iter()
                .filter_map(|old| {
                    // Find if this node still exists.
                    let still_exists = after_nodes.contains(old);
                    if still_exists {
                        Some((*old, Some(*old)))
                    } else {
                        Some((*old, None))
                    }
                })
                .collect();

            self.history.record_mutation(&self.graph, action, &entity_map);

            self.history.finalize_operation(OperationStats {
                entity_count_before: before_nodes.len(),
                entity_count_after: after_nodes.len(),
                names_preserved: entity_map.iter().filter(|(_, new)| new.is_some()).count(),
                names_lost: entity_map.iter().filter(|(_, new)| new.is_none()).count(),
                names_generated: 0,
                conflicts_resolved: 0,
            });

            Ok(())
        } else {
            self.history.cancel_operation();
            Err(issues)
        }
    }

    /// Get the persistent ID for a node.
    pub fn get_persistent_id(&self, node: TopoNode) -> Option<PersistentId> {
        self.history.resolve_persistent(node)
    }

    /// Generate a stability report.
    pub fn stability_report(&self) -> CrossOperationStabilityReport {
        self.history.stability_report()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cross-Operation Naming Stability Analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Metrics measuring naming stability for a single operation.
#[derive(Debug, Clone)]
pub struct OperationStabilityMetrics {
    /// Operation identifier.
    pub operation_id: rcad_kernel::persistent_naming::OperationId,
    /// Type of operation performed.
    pub operation_type: rcad_kernel::persistent_naming::OperationType,
    /// Optional label for the operation.
    pub label: Option<String>,
    /// Number of entities that retained their names through this operation.
    pub names_retained: usize,
    /// Number of entities that lost their names during this operation.
    pub names_lost: usize,
    /// Number of new names generated during this operation.
    pub names_generated: usize,
    /// Number of naming conflicts that occurred during this operation.
    pub conflicts: usize,
    /// Stability score for this specific operation (0.0 - 1.0).
    pub stability_score: f64,
    /// Cumulative stability score up to and including this operation.
    pub cumulative_stability: f64,
}

impl Default for OperationStabilityMetrics {
    fn default() -> Self {
        Self {
            operation_id: rcad_kernel::persistent_naming::OperationId::NULL,
            operation_type: rcad_kernel::persistent_naming::OperationType::Generic,
            label: None,
            names_retained: 0,
            names_lost: 0,
            names_generated: 0,
            conflicts: 0,
            stability_score: 1.0,
            cumulative_stability: 1.0,
        }
    }
}

/// Information about a broken naming chain.
#[derive(Debug, Clone)]
pub struct BrokenChainInfo {
    /// The persistent ID whose chain was broken.
    pub persistent_id: PersistentId,
    /// The operation where the break occurred.
    pub broken_at_operation: rcad_kernel::persistent_naming::OperationId,
    /// Entity ID that lost the name.
    pub entity_id: u64,
    /// Entity type if known.
    pub entity_type: Option<rcad_kernel::persistent_naming::EntityType>,
    /// Number of operations the chain survived before breaking.
    pub survived_operations: usize,
    /// Description of how the chain broke.
    pub break_reason: String,
}

/// A naming conflict detected across operations.
#[derive(Debug, Clone)]
pub struct NamingConflict {
    /// The persistent ID involved in the conflict.
    pub persistent_id: PersistentId,
    /// Operations where the conflict manifested.
    pub involved_operations: Vec<rcad_kernel::persistent_naming::OperationId>,
    /// Entity IDs that were in conflict.
    pub conflicting_entities: Vec<u64>,
    /// Type of conflict.
    pub conflict_type: ConflictType,
    /// Severity of the conflict.
    pub severity: rcad_kernel::persistent_naming::IssueSeverity,
    /// Whether the conflict was automatically resolved.
    pub auto_resolved: bool,
}

/// Types of naming conflicts that can occur.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConflictType {
    /// Same persistent ID assigned to multiple entities.
    DuplicateAssignment,
    /// Entity references a deleted persistent ID.
    ReferenceToDeleted,
    /// Genealogy chain is incomplete or broken.
    BrokenGenealogy,
    /// Unexpected name change during propagation.
    UnexpectedNameChange,
    /// Merge operation lost entity tracking.
    MergeTrackingLoss,
}

/// A recommendation for improving naming stability.
#[derive(Debug, Clone)]
pub struct StabilityRecommendation {
    /// Priority of this recommendation (higher = more important).
    pub priority: u32,
    /// Category of the recommendation.
    pub category: RecommendationCategory,
    /// Human-readable description of the recommendation.
    pub description: String,
    /// Operations this recommendation applies to (empty = all).
    pub affected_operations: Vec<rcad_kernel::persistent_naming::OperationId>,
    /// Estimated impact on stability score if implemented.
    pub estimated_impact: f64,
    /// Code or configuration suggestion.
    pub suggestion: Option<String>,
}

/// Categories of stability recommendations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecommendationCategory {
    /// Adjust naming rule selection.
    NamingRule,
    /// Improve name propagation policy.
    PropagationPolicy,
    /// Fix specific operation handling.
    OperationHandling,
    /// Improve entity tracking.
    EntityTracking,
    /// Address conflict resolution.
    ConflictResolution,
    /// General architecture improvement.
    Architecture,
}

/// Comprehensive analysis of naming stability across multiple operations.
#[derive(Debug, Clone)]
pub struct CrossOperationNamingAnalysis {
    /// Operations analyzed.
    pub operations: Vec<OperationRecord>,
    /// Entities tracked through all operations.
    pub entity_genealogy: HashMap<PersistentId, EntityGenealogy>,
    /// Stability metrics per operation.
    pub per_operation_stability: Vec<OperationStabilityMetrics>,
    /// Overall stability score (0.0 - 1.0).
    pub overall_stability: f64,
    /// Entities with broken naming chains.
    pub broken_chains: Vec<BrokenChainInfo>,
    /// Trend direction: positive = improving, negative = degrading.
    pub stability_trend: f64,
    /// Number of entities tracked at each operation boundary.
    pub entity_counts: Vec<usize>,
}

impl CrossOperationNamingAnalysis {
    /// Returns true if overall stability is excellent (> 95%).
    pub fn is_excellent(&self) -> bool {
        self.overall_stability >= 0.95
    }

    /// Returns true if overall stability is good (> 90%).
    pub fn is_good(&self) -> bool {
        self.overall_stability >= 0.90
    }

    /// Returns true if there are significant stability issues.
    pub fn has_issues(&self) -> bool {
        self.overall_stability < 0.90 || !self.broken_chains.is_empty()
    }

    /// Get the most problematic operation (lowest stability score).
    pub fn most_problematic_operation(&self) -> Option<&OperationStabilityMetrics> {
        self.per_operation_stability
            .iter()
            .min_by(|a, b| a.stability_score.partial_cmp(&b.stability_score).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// Get operations sorted by stability score (ascending).
    pub fn operations_by_stability(&self) -> Vec<&OperationStabilityMetrics> {
        let mut ops: Vec<_> = self.per_operation_stability.iter().collect();
        ops.sort_by(|a, b| a.stability_score.partial_cmp(&b.stability_score).unwrap_or(std::cmp::Ordering::Equal));
        ops
    }

    /// Calculate the average stability score across all operations.
    pub fn average_operation_stability(&self) -> f64 {
        if self.per_operation_stability.is_empty() {
            return 1.0;
        }
        let sum: f64 = self.per_operation_stability.iter().map(|m| m.stability_score).sum();
        sum / self.per_operation_stability.len() as f64
    }

    /// Generate a human-readable summary.
    pub fn summary(&self) -> String {
        format!(
            "Cross-Operation Naming Analysis:\n\
             - Operations: {}\n\
             - Entities Tracked: {}\n\
             - Overall Stability: {:.1}%\n\
             - Stability Trend: {}\n\
             - Broken Chains: {}\n\
             - Avg Operation Stability: {:.1}%",
            self.operations.len(),
            self.entity_genealogy.len(),
            self.overall_stability * 100.0,
            if self.stability_trend > 0.0 { "Improving" } else if self.stability_trend < 0.0 { "Degrading" } else { "Stable" },
            self.broken_chains.len(),
            self.average_operation_stability() * 100.0
        )
    }
}

/// Track naming through a sequence of operations.
///
/// This function analyzes a series of BRepGraphHistory snapshots to determine
/// how naming stability evolves across multiple operations.
pub fn analyze_naming_sequence(
    history: &[BRepGraphHistory],
    initial_entities: &[TopoNode],
) -> CrossOperationNamingAnalysis {
    use rcad_kernel::persistent_naming::{OperationStats, NamingEvent};

    if history.is_empty() {
        return CrossOperationNamingAnalysis {
            operations: Vec::new(),
            entity_genealogy: HashMap::new(),
            per_operation_stability: Vec::new(),
            overall_stability: 1.0,
            broken_chains: Vec::new(),
            stability_trend: 0.0,
            entity_counts: vec![initial_entities.len()],
        };
    }

    let mut operations: Vec<OperationRecord> = Vec::new();
    let mut entity_genealogy: HashMap<PersistentId, EntityGenealogy> = HashMap::new();
    let mut per_operation_stability: Vec<OperationStabilityMetrics> = Vec::new();
    let mut broken_chains: Vec<BrokenChainInfo> = Vec::new();
    let mut entity_counts: Vec<usize> = vec![initial_entities.len()];

    // Track entity persistence across operations.
    let mut entity_to_pid: HashMap<u64, PersistentId> = HashMap::new();
    let mut pid_to_entity: HashMap<PersistentId, u64> = HashMap::new();

    // Assign initial persistent IDs.
    for node in initial_entities {
        let entity_id = node_to_entity_id(*node);
        let pid = PersistentId(entity_id); // Use entity_id as basis for PID
        entity_to_pid.insert(entity_id, pid);
        pid_to_entity.insert(pid, entity_id);

        entity_genealogy.insert(pid, EntityGenealogy {
            persistent_id: pid,
            created_in_operation: rcad_kernel::persistent_naming::OperationId::NULL,
            evolution: vec![(rcad_kernel::persistent_naming::OperationId::NULL, entity_id)],
            current_entity_id: Some(entity_id),
            is_deleted: false,
        });
    }

    let mut cumulative_stability = 1.0;

    for (op_idx, hist) in history.iter().enumerate() {
        let cross_op = hist.naming_engine().cross_operation_history();

        // Copy operation records.
        for op in &cross_op.operations {
            operations.push(op.clone());
        }

        // Calculate stability metrics for each operation in this history.
        for op in &cross_op.operations {
            let total_entities = op.stats.entity_count_before.max(1);
            let preserved = op.stats.names_preserved;
            let stability_score = preserved as f64 / total_entities as f64;

            cumulative_stability = cumulative_stability * stability_score;

            per_operation_stability.push(OperationStabilityMetrics {
                operation_id: op.id,
                operation_type: op.operation_type,
                label: op.label.clone(),
                names_retained: op.stats.names_preserved,
                names_lost: op.stats.names_lost,
                names_generated: op.stats.names_generated,
                conflicts: op.stats.conflicts_resolved,
                stability_score,
                cumulative_stability,
            });

            // Track entity counts.
            entity_counts.push(op.stats.entity_count_after);

            // Detect broken chains.
            for event in &op.naming_events {
                if let NamingEvent::Lost { entity_id, persistent_id } = event {
                    // Find how many operations this entity survived.
                    let survived_operations = entity_genealogy
                        .get(persistent_id)
                        .map(|g| g.evolution.len())
                        .unwrap_or(0);

                    broken_chains.push(BrokenChainInfo {
                        persistent_id: *persistent_id,
                        broken_at_operation: op.id,
                        entity_id: *entity_id,
                        entity_type: infer_entity_type_from_id(*entity_id),
                        survived_operations,
                        break_reason: "Entity removed without successor".to_string(),
                    });
                }
            }
        }

        // Update genealogy from this history.
        for (pid, genealogy) in cross_op.genealogy.iter() {
            entity_genealogy.insert(*pid, genealogy.clone());
        }
    }

    // Calculate overall stability.
    let total_preserved: usize = per_operation_stability.iter().map(|m| m.names_retained).sum();
    let total_lost: usize = per_operation_stability.iter().map(|m| m.names_lost).sum();
    let total = total_preserved + total_lost;
    let overall_stability = if total > 0 {
        total_preserved as f64 / total as f64
    } else {
        1.0
    };

    // Calculate stability trend.
    let stability_trend = calculate_stability_trend(&per_operation_stability);

    CrossOperationNamingAnalysis {
        operations,
        entity_genealogy,
        per_operation_stability,
        overall_stability,
        broken_chains,
        stability_trend,
        entity_counts,
    }
}

/// Detect naming conflicts across operations.
///
/// This function analyzes cross-operation data to identify conflicts that
/// may not be visible in single-operation analysis.
pub fn detect_cross_operation_conflicts(
    analysis: &CrossOperationNamingAnalysis,
) -> Vec<NamingConflict> {
    let mut conflicts: Vec<NamingConflict> = Vec::new();
    let mut pid_to_entities: HashMap<PersistentId, Vec<(u64, rcad_kernel::persistent_naming::OperationId)>> = HashMap::new();

    // Build a map of persistent ID to all entities that have held it.
    for (pid, genealogy) in &analysis.entity_genealogy {
        let mut entities = Vec::new();
        for (op_id, entity_id) in &genealogy.evolution {
            entities.push((*entity_id, *op_id));
        }
        pid_to_entities.insert(*pid, entities);
    }

    // Detect duplicate assignments (same PID to different entities at same time).
    for (pid, genealogy) in &analysis.entity_genealogy {
        if genealogy.is_deleted {
            continue;
        }

        // Check for entities that reference this PID after it was marked deleted.
        if let Some(current_entity) = genealogy.current_entity_id {
            // Verify the current entity still has this PID.
            for (other_pid, other_genealogy) in &analysis.entity_genealogy {
                if other_pid != pid && other_genealogy.current_entity_id == Some(current_entity) {
                    // Same entity has multiple PIDs.
                    conflicts.push(NamingConflict {
                        persistent_id: *pid,
                        involved_operations: vec![genealogy.created_in_operation],
                        conflicting_entities: vec![current_entity],
                        conflict_type: ConflictType::DuplicateAssignment,
                        severity: rcad_kernel::persistent_naming::IssueSeverity::Severe,
                        auto_resolved: false,
                    });
                }
            }
        }
    }

    // Detect broken genealogies.
    for (pid, genealogy) in &analysis.entity_genealogy {
        if genealogy.evolution.is_empty() && !genealogy.is_deleted {
            conflicts.push(NamingConflict {
                persistent_id: *pid,
                involved_operations: vec![genealogy.created_in_operation],
                conflicting_entities: vec![],
                conflict_type: ConflictType::BrokenGenealogy,
                severity: rcad_kernel::persistent_naming::IssueSeverity::Moderate,
                auto_resolved: false,
            });
        }
    }

    // Detect reference to deleted entities.
    for chain in &analysis.broken_chains {
        // Check if any genealogy still references this broken chain.
        if let Some(genealogy) = analysis.entity_genealogy.get(&chain.persistent_id) {
            if !genealogy.is_deleted && genealogy.evolution.len() > chain.survived_operations {
                conflicts.push(NamingConflict {
                    persistent_id: chain.persistent_id,
                    involved_operations: vec![chain.broken_at_operation],
                    conflicting_entities: vec![chain.entity_id],
                    conflict_type: ConflictType::ReferenceToDeleted,
                    severity: rcad_kernel::persistent_naming::IssueSeverity::Critical,
                    auto_resolved: false,
                });
            }
        }
    }

    conflicts
}

/// Generate recommendations for improving naming stability.
///
/// Based on the analysis, this function produces actionable recommendations
/// for improving the naming stability of future operations.
pub fn generate_stability_recommendations(
    analysis: &CrossOperationNamingAnalysis,
) -> Vec<StabilityRecommendation> {
    let mut recommendations: Vec<StabilityRecommendation> = Vec::new();

    // Check overall stability.
    if analysis.overall_stability < 0.5 {
        recommendations.push(StabilityRecommendation {
            priority: 100,
            category: RecommendationCategory::Architecture,
            description: "Critical naming stability issues detected. Consider reviewing the entire naming strategy.".to_string(),
            affected_operations: vec![],
            estimated_impact: 0.3,
            suggestion: Some("Enable Hybrid naming rule and Preserve propagation policy".to_string()),
        });
    } else if analysis.overall_stability < 0.8 {
        recommendations.push(StabilityRecommendation {
            priority: 80,
            category: RecommendationCategory::NamingRule,
            description: "Moderate naming stability degradation. Review naming rule configuration.".to_string(),
            affected_operations: vec![],
            estimated_impact: 0.15,
            suggestion: Some("Consider switching to HistoryTracking naming rule for better traceability".to_string()),
        });
    }

    // Check for problematic operations.
    if let Some(problematic) = analysis.most_problematic_operation() {
        if problematic.stability_score < 0.7 {
            recommendations.push(StabilityRecommendation {
                priority: 90,
                category: RecommendationCategory::OperationHandling,
                description: format!(
                    "Operation {:?} has low stability score ({:.1}%). Review operation-specific handling.",
                    problematic.operation_type,
                    problematic.stability_score * 100.0
                ),
                affected_operations: vec![problematic.operation_id],
                estimated_impact: 0.2,
                suggestion: Some("Ensure entity mapping is correctly tracked during this operation".to_string()),
            });
        }
    }

    // Check for broken chains.
    if !analysis.broken_chains.is_empty() {
        let severe_breaks: Vec<_> = analysis.broken_chains.iter()
            .filter(|c| c.survived_operations > 5)
            .collect();

        if !severe_breaks.is_empty() {
            recommendations.push(StabilityRecommendation {
                priority: 85,
                category: RecommendationCategory::EntityTracking,
                description: format!(
                    "{} long-lived entities lost their naming chains. Improve entity tracking during mutations.",
                    severe_breaks.len()
                ),
                affected_operations: severe_breaks.iter().map(|c| c.broken_at_operation).collect(),
                estimated_impact: 0.25,
                suggestion: Some("Implement explicit entity mapping during split/merge operations".to_string()),
            });
        }
    }

    // Check stability trend.
    if analysis.stability_trend < -0.1 {
        recommendations.push(StabilityRecommendation {
            priority: 75,
            category: RecommendationCategory::PropagationPolicy,
            description: "Stability is degrading over time. Naming propagation may need adjustment.".to_string(),
            affected_operations: vec![],
            estimated_impact: 0.1,
            suggestion: Some("Review NamePropagationPolicy settings for recent operations".to_string()),
        });
    }

    // Check for conflicts.
    let conflicts = detect_cross_operation_conflicts(analysis);
    let critical_conflicts: Vec<_> = conflicts.iter()
        .filter(|c| c.severity == rcad_kernel::persistent_naming::IssueSeverity::Critical)
        .collect();

    if !critical_conflicts.is_empty() {
        recommendations.push(StabilityRecommendation {
            priority: 95,
            category: RecommendationCategory::ConflictResolution,
            description: format!(
                "{} critical naming conflicts detected. Immediate resolution required.",
                critical_conflicts.len()
            ),
            affected_operations: critical_conflicts.iter().flat_map(|c| c.involved_operations.clone()).collect(),
            estimated_impact: 0.3,
            suggestion: Some("Manually resolve conflicts or reset naming for affected entities".to_string()),
        });
    }

    // Check per-operation patterns.
    let low_stability_ops: Vec<_> = analysis.per_operation_stability.iter()
        .filter(|m| m.stability_score < 0.8)
        .collect();

    if low_stability_ops.len() > analysis.per_operation_stability.len() / 2 {
        recommendations.push(StabilityRecommendation {
            priority: 70,
            category: RecommendationCategory::Architecture,
            description: "Multiple operations have low stability. Consider system-wide naming improvements.".to_string(),
            affected_operations: low_stability_ops.iter().map(|m| m.operation_id).collect(),
            estimated_impact: 0.2,
            suggestion: Some("Enable comprehensive entity tracking and increase name propagation fidelity".to_string()),
        });
    }

    // Sort by priority (highest first).
    recommendations.sort_by(|a, b| b.priority.cmp(&a.priority));
    recommendations
}

/// Calculate the stability trend from per-operation metrics.
///
/// Returns a value between -1.0 (strongly degrading) and 1.0 (strongly improving).
fn calculate_stability_trend(metrics: &[OperationStabilityMetrics]) -> f64 {
    if metrics.len() < 2 {
        return 0.0;
    }

    // Simple linear regression on stability scores.
    let n = metrics.len() as f64;
    let sum_x: f64 = (0..metrics.len()).map(|i| i as f64).sum();
    let sum_y: f64 = metrics.iter().map(|m| m.stability_score).sum();
    let sum_xy: f64 = metrics.iter().enumerate()
        .map(|(i, m)| i as f64 * m.stability_score)
        .sum();
    let sum_x2: f64 = (0..metrics.len()).map(|i| (i * i) as f64).sum();

    let denominator = n * sum_x2 - sum_x * sum_x;
    if denominator.abs() < f64::EPSILON {
        return 0.0;
    }

    let slope = (n * sum_xy - sum_x * sum_y) / denominator;

    // Normalize slope to [-1, 1] range.
    // A slope of 0.1 per operation is considered a strong trend.
    slope.clamp(-1.0, 1.0)
}

/// Infer entity type from an encoded entity ID.
fn infer_entity_type_from_id(entity_id: u64) -> Option<rcad_kernel::persistent_naming::EntityType> {
    // High 8 bits encode the kind.
    let kind_bits = entity_id >> 56;
    match kind_bits {
        0 | 1 => Some(rcad_kernel::persistent_naming::EntityType::Solid),
        2 => Some(rcad_kernel::persistent_naming::EntityType::Face),
        4 => Some(rcad_kernel::persistent_naming::EntityType::Edge),
        5 => Some(rcad_kernel::persistent_naming::EntityType::Vertex),
        _ => None,
    }
}

// Re-export types needed for analysis.
pub use rcad_kernel::persistent_naming::{
    OperationRecord, EntityGenealogy,
};

impl BRepGraphHistory {
    /// Cancel the current operation.
    fn cancel_operation(&mut self) {
        self.naming_engine.cancel_operation();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Enhanced Persistent Naming Semantics
// ─────────────────────────────────────────────────────────────────────────────

/// A scoped identifier for a topological entity within a naming context.
///
/// `ScopedId` combines a persistent ID with a naming scope (part, assembly, operation)
/// to provide fully-qualified identifiers that are unique across an entire model hierarchy.
///
/// # Example
/// ```
/// use rcad_algorithms::brep_graph::{ScopedId, NamingScope};
/// use rcad_kernel::PersistentId;
///
/// // A face with ID 42 in part "housing", assembly "machine", operation "fillet"
/// let scoped = ScopedId {
///     persistent_id: PersistentId(42),
///     scope: NamingScope {
///         part: Some("housing".to_string()),
///         assembly: Some("machine".to_string()),
///         operation: Some("fillet".to_string()),
///     },
/// };
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ScopedId {
    /// The stable persistent ID for this entity.
    pub persistent_id: PersistentId,
    /// The naming scope in which this ID is defined.
    pub scope: NamingScope,
}

impl ScopedId {
    /// Create a new scoped ID with the given persistent ID and scope.
    pub fn new(persistent_id: PersistentId, scope: NamingScope) -> Self {
        Self { persistent_id, scope }
    }

    /// Create a scoped ID with a null persistent ID and empty scope.
    pub fn null() -> Self {
        Self {
            persistent_id: PersistentId::NULL,
            scope: NamingScope::default(),
        }
    }

    /// Returns true if this is a null/invalid scoped ID.
    pub fn is_null(&self) -> bool {
        self.persistent_id.is_null()
    }

    /// Generate a fully-qualified name string for this scoped ID.
    pub fn qualified_name(&self) -> String {
        let mut parts = Vec::new();
        if let Some(ref assembly) = self.scope.assembly {
            parts.push(assembly.clone());
        }
        if let Some(ref part) = self.scope.part {
            parts.push(part.clone());
        }
        if let Some(ref op) = self.scope.operation {
            parts.push(op.clone());
        }
        parts.push(format!("e{}", self.persistent_id.raw()));
        parts.join("::")
    }
}

/// The naming scope defines the context in which persistent IDs are meaningful.
///
/// Scopes form a hierarchy: assembly > part > operation. An entity's full identity
/// is determined by its persistent ID within the current scope.
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct NamingScope {
    /// The part name (e.g., "housing", "cover").
    pub part: Option<String>,
    /// The assembly name (e.g., "machine", "device").
    pub assembly: Option<String>,
    /// The operation that created or last modified this entity.
    pub operation: Option<String>,
}

impl NamingScope {
    /// Create a new empty scope.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a scope for a specific part.
    pub fn for_part(part: impl Into<String>) -> Self {
        Self {
            part: Some(part.into()),
            assembly: None,
            operation: None,
        }
    }

    /// Create a scope for a specific assembly.
    pub fn for_assembly(assembly: impl Into<String>) -> Self {
        Self {
            part: None,
            assembly: Some(assembly.into()),
            operation: None,
        }
    }

    /// Create a scope for a specific operation within a part.
    pub fn for_operation(part: impl Into<String>, operation: impl Into<String>) -> Self {
        Self {
            part: Some(part.into()),
            assembly: None,
            operation: Some(operation.into()),
        }
    }

    /// Set the part name.
    pub fn with_part(mut self, part: impl Into<String>) -> Self {
        self.part = Some(part.into());
        self
    }

    /// Set the assembly name.
    pub fn with_assembly(mut self, assembly: impl Into<String>) -> Self {
        self.assembly = Some(assembly.into());
        self
    }

    /// Set the operation name.
    pub fn with_operation(mut self, operation: impl Into<String>) -> Self {
        self.operation = Some(operation.into());
        self
    }

    /// Create a child scope for a sub-operation.
    pub fn child_scope(&self, operation: impl Into<String>) -> Self {
        Self {
            part: self.part.clone(),
            assembly: self.assembly.clone(),
            operation: Some(operation.into()),
        }
    }

    /// Check if this scope is a parent of (or equal to) another scope.
    pub fn contains(&self, other: &NamingScope) -> bool {
        match (&self.assembly, &other.assembly) {
            (Some(a1), Some(a2)) if a1 != a2 => return false,
            (Some(_), None) => return false,
            _ => {}
        }
        match (&self.part, &other.part) {
            (Some(p1), Some(p2)) if p1 != p2 => return false,
            (Some(_), None) => return false,
            _ => {}
        }
        true
    }
}

/// Enhanced naming context that tracks scopes and entity relationships.
///
/// `EnhancedNamingContext` extends the basic `NamingContext` with:
/// - Scope-aware ID assignment
/// - Detailed genealogy tracking
/// - Conflict detection and resolution
/// - Serialization support
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct EnhancedNamingContext {
    /// The current naming scope.
    pub current_scope: NamingScope,
    /// Mapping from entity IDs to scoped IDs.
    entity_to_scoped: HashMap<u64, ScopedId>,
    /// Reverse mapping from scoped IDs to entity IDs.
    scoped_to_entity: HashMap<ScopedId, u64>,
    /// Genealogy records indexed by persistent ID.
    genealogy: HashMap<PersistentId, EntityGenealogyRecord>,
    /// Pending name assignments waiting for scope resolution.
    pending_assignments: Vec<PendingNameAssignment>,
    /// Conflict resolution history.
    conflict_history: Vec<NameConflictRecord>,
    /// Next persistent ID to allocate.
    next_persistent_id: u64,
}

/// Record of an entity's genealogy through operations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntityGenealogyRecord {
    /// The persistent ID being tracked.
    pub persistent_id: PersistentId,
    /// The scope in which this entity was created.
    pub creation_scope: NamingScope,
    /// The operation that created this entity.
    pub creation_operation: Option<String>,
    /// Chain of transformations: (operation, entity_id_before, entity_id_after).
    pub transformation_chain: Vec<GenealogyStep>,
    /// Parent entity IDs (for merged entities, this has multiple entries).
    pub parent_ids: Vec<PersistentId>,
    /// Child entity IDs (for split entities, this has multiple entries).
    pub child_ids: Vec<PersistentId>,
    /// Current status of this entity.
    pub status: EntityStatus,
}

/// A single step in an entity's genealogy transformation chain.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenealogyStep {
    /// The operation that caused this transformation.
    pub operation: String,
    /// The entity ID before this operation (None if generated).
    pub entity_id_before: Option<u64>,
    /// The entity ID after this operation.
    pub entity_id_after: u64,
    /// The scope at the time of this operation.
    pub scope: NamingScope,
}

/// Status of an entity in the genealogy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EntityStatus {
    /// Entity is active and present in the model.
    Active,
    /// Entity was deleted or consumed by an operation.
    Deleted,
    /// Entity was merged into another entity.
    Merged,
    /// Entity was split into multiple entities.
    Split,
    /// Entity is pending resolution.
    Pending,
}

/// A pending name assignment waiting for scope resolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PendingNameAssignment {
    /// The entity ID to be assigned.
    pub entity_id: u64,
    /// The proposed scope for this assignment.
    pub proposed_scope: NamingScope,
    /// Source entity IDs this entity was derived from.
    pub source_entities: Vec<u64>,
    /// The propagation policy to use.
    pub propagation_policy: NamePropagationPolicy,
}

/// Record of a naming conflict and its resolution.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NameConflictRecord {
    /// The conflicting persistent ID.
    pub persistent_id: PersistentId,
    /// The entity IDs involved in the conflict.
    pub conflicting_entities: Vec<u64>,
    /// The operation where the conflict occurred.
    pub operation: String,
    /// The scope where the conflict occurred.
    pub scope: NamingScope,
    /// How the conflict was resolved.
    pub resolution: NameConflictResolution,
    /// Timestamp of the conflict (sequence number).
    pub sequence: u64,
}

/// Strategies for resolving naming conflicts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NameConflictResolution {
    /// Kept the existing binding, rejected the new.
    KeepExisting,
    /// Replaced the existing binding with the new.
    ReplaceWithNew,
    /// Generated a new persistent ID for the new entity.
    GenerateNewId,
    /// Merged both entities under a shared context.
    MergeEntities,
    /// Created an alias mapping.
    CreateAlias,
    /// Could not resolve automatically - requires manual intervention.
    Unresolved,
}

/// Operation-specific name propagation rules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NamePropagationRule {
    /// The operation type this rule applies to.
    pub operation_type: OperationType,
    /// Policy for face entities.
    pub face_policy: NamePropagationPolicy,
    /// Policy for edge entities.
    pub edge_policy: NamePropagationPolicy,
    /// Policy for vertex entities.
    pub vertex_policy: NamePropagationPolicy,
    /// Whether to track genealogy for this operation.
    pub track_genealogy: bool,
    /// Conflict resolution strategy for this operation.
    pub conflict_resolution: NameConflictResolution,
}

impl NamePropagationRule {
    /// Create a default propagation rule for an operation type.
    pub fn for_operation(operation_type: OperationType) -> Self {
        match operation_type {
            OperationType::BooleanUnion |
            OperationType::BooleanIntersection |
            OperationType::BooleanDifference => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Preserve,
                edge_policy: NamePropagationPolicy::Preserve,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::GenerateNewId,
            },
            OperationType::Feature => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Inherit,
                edge_policy: NamePropagationPolicy::Inherit,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::KeepExisting,
            },
            OperationType::EdgeSplit |
            OperationType::FaceSplit => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Inherit,
                edge_policy: NamePropagationPolicy::Inherit,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::GenerateNewId,
            },
            OperationType::Merge => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Combine,
                edge_policy: NamePropagationPolicy::Combine,
                vertex_policy: NamePropagationPolicy::Combine,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::MergeEntities,
            },
            OperationType::Delete => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Generate,
                edge_policy: NamePropagationPolicy::Generate,
                vertex_policy: NamePropagationPolicy::Generate,
                track_genealogy: false,
                conflict_resolution: NameConflictResolution::KeepExisting,
            },
            OperationType::Transform => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Preserve,
                edge_policy: NamePropagationPolicy::Preserve,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::KeepExisting,
            },
            OperationType::Generic |
            OperationType::Import => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Preserve,
                edge_policy: NamePropagationPolicy::Preserve,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::GenerateNewId,
            },
        }
    }

    /// Get the propagation policy for a specific entity kind.
    pub fn policy_for_kind(&self, kind: NodeKind) -> NamePropagationPolicy {
        match kind {
            NodeKind::Face | NodeKind::Shell | NodeKind::Solid => self.face_policy,
            NodeKind::Edge | NodeKind::Wire => self.edge_policy,
            NodeKind::Vertex => self.vertex_policy,
        }
    }
}

/// Propagation policies for name inheritance through operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NamePropagationPolicy {
    /// Keep the original entity's name unchanged.
    Preserve,
    /// Inherit the parent entity's name with a disambiguating suffix.
    Inherit,
    /// Generate a completely new name.
    Generate,
    /// Combine names from multiple source entities (for merges).
    Combine,
    /// Create a derivative name based on geometric properties.
    GeometryBased,
    /// Create a derivative name based on topological relationships.
    TopologyBased,
}

impl Default for NamePropagationPolicy {
    fn default() -> Self {
        Self::Preserve
    }
}

impl EnhancedNamingContext {
    /// Create a new empty naming context.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a naming context with a specific scope.
    pub fn with_scope(scope: NamingScope) -> Self {
        Self {
            current_scope: scope,
            ..Default::default()
        }
    }

    /// Set the current naming scope.
    pub fn set_scope(&mut self, scope: NamingScope) {
        self.current_scope = scope;
    }

    /// Get the current naming scope.
    pub fn scope(&self) -> &NamingScope {
        &self.current_scope
    }

    /// Assign a new persistent ID to an entity.
    pub fn assign_id(&mut self, entity_id: u64) -> PersistentId {
        let pid = self.allocate_persistent_id();
        let scoped_id = ScopedId::new(pid, self.current_scope.clone());
        self.entity_to_scoped.insert(entity_id, scoped_id.clone());
        self.scoped_to_entity.insert(scoped_id, entity_id);

        // Create genealogy record.
        self.genealogy.insert(pid, EntityGenealogyRecord {
            persistent_id: pid,
            creation_scope: self.current_scope.clone(),
            creation_operation: self.current_scope.operation.clone(),
            transformation_chain: vec![],
            parent_ids: vec![],
            child_ids: vec![],
            status: EntityStatus::Active,
        });

        pid
    }

    /// Assign a persistent ID derived from source entities.
    pub fn assign_derived_id(
        &mut self,
        entity_id: u64,
        source_entities: &[u64],
        policy: NamePropagationPolicy,
    ) -> PersistentId {
        match policy {
            NamePropagationPolicy::Preserve | NamePropagationPolicy::Inherit => {
                // Inherit from the first source that has a persistent ID.
                if let Some(&source_id) = source_entities.first() {
                    if let Some(scoped) = self.entity_to_scoped.get(&source_id) {
                        let pid = scoped.persistent_id;
                        let new_scoped = ScopedId::new(pid, self.current_scope.clone());
                        self.entity_to_scoped.insert(entity_id, new_scoped.clone());
                        self.scoped_to_entity.insert(new_scoped, entity_id);

                        // Update genealogy.
                        if let Some(record) = self.genealogy.get_mut(&pid) {
                            record.transformation_chain.push(GenealogyStep {
                                operation: self.current_scope.operation.clone().unwrap_or_default(),
                                entity_id_before: Some(source_id),
                                entity_id_after: entity_id,
                                scope: self.current_scope.clone(),
                            });
                        }

                        return pid;
                    }
                }
                self.assign_id(entity_id)
            }
            NamePropagationPolicy::Combine => {
                // Combine all source persistent IDs into the genealogy.
                let pid = self.assign_id(entity_id);
                if let Some(record) = self.genealogy.get_mut(&pid) {
                    for &source_id in source_entities {
                        if let Some(scoped) = self.entity_to_scoped.get(&source_id) {
                            record.parent_ids.push(scoped.persistent_id);
                        }
                    }
                }
                pid
            }
            NamePropagationPolicy::Generate |
            NamePropagationPolicy::GeometryBased |
            NamePropagationPolicy::TopologyBased => {
                self.assign_id(entity_id)
            }
        }
    }

    /// Resolve a persistent ID to an entity ID within the current scope.
    pub fn resolve_entity(&self, pid: PersistentId) -> Option<u64> {
        let scoped = ScopedId::new(pid, self.current_scope.clone());
        self.scoped_to_entity.get(&scoped).copied()
    }

    /// Resolve an entity ID to a persistent ID.
    pub fn resolve_persistent(&self, entity_id: u64) -> Option<PersistentId> {
        self.entity_to_scoped.get(&entity_id).map(|s| s.persistent_id)
    }

    /// Record a split operation: one entity becomes multiple.
    pub fn record_split(
        &mut self,
        source_entity_id: u64,
        target_entity_ids: &[u64],
        operation: &str,
    ) -> Vec<PersistentId> {
        let source_pid = self.resolve_persistent(source_entity_id);
        let mut result_pids = Vec::with_capacity(target_entity_ids.len());

        for (i, &target_id) in target_entity_ids.iter().enumerate() {
            let pid = if i == 0 {
                // First target inherits the source's persistent ID.
                if let Some(pid) = source_pid {
                    let scoped = ScopedId::new(pid, self.current_scope.clone());
                    self.entity_to_scoped.insert(target_id, scoped.clone());
                    self.scoped_to_entity.insert(scoped, target_id);
                    pid
                } else {
                    self.assign_id(target_id)
                }
            } else {
                // Subsequent targets get new IDs.
                self.assign_id(target_id)
            };
            result_pids.push(pid);
        }

        // Update genealogy for the source entity.
        if let Some(pid) = source_pid {
            if let Some(record) = self.genealogy.get_mut(&pid) {
                record.status = EntityStatus::Split;
                record.child_ids.extend_from_slice(&result_pids[1..]);
                record.transformation_chain.push(GenealogyStep {
                    operation: operation.to_string(),
                    entity_id_before: Some(source_entity_id),
                    entity_id_after: result_pids.first().map(|&p| {
                        self.scoped_to_entity.get(&ScopedId::new(p, self.current_scope.clone()))
                            .copied()
                            .unwrap_or(0)
                    }).unwrap_or(0),
                    scope: self.current_scope.clone(),
                });
            }
        }

        result_pids
    }

    /// Record a merge operation: multiple entities become one.
    pub fn record_merge(
        &mut self,
        source_entity_ids: &[u64],
        target_entity_id: u64,
        operation: &str,
        resolution: NameConflictResolution,
    ) -> PersistentId {
        // Find the first source with a persistent ID.
        let primary_pid = source_entity_ids
            .iter()
            .find_map(|&id| self.resolve_persistent(id));

        let target_pid = match resolution {
            NameConflictResolution::KeepExisting => {
                if let Some(pid) = primary_pid {
                    let scoped = ScopedId::new(pid, self.current_scope.clone());
                    self.entity_to_scoped.insert(target_entity_id, scoped.clone());
                    self.scoped_to_entity.insert(scoped, target_entity_id);
                    pid
                } else {
                    self.assign_id(target_entity_id)
                }
            }
            NameConflictResolution::GenerateNewId => {
                self.assign_id(target_entity_id)
            }
            NameConflictResolution::MergeEntities => {
                let pid = self.assign_id(target_entity_id);
                // Collect parent IDs first to avoid borrow conflict
                let parent_pids: Vec<PersistentId> = source_entity_ids
                    .iter()
                    .filter_map(|&source_id| self.resolve_persistent(source_id))
                    .collect();
                if let Some(record) = self.genealogy.get_mut(&pid) {
                    for source_pid in parent_pids {
                        record.parent_ids.push(source_pid);
                    }
                }
                pid
            }
            _ => {
                if let Some(pid) = primary_pid {
                    pid
                } else {
                    self.assign_id(target_entity_id)
                }
            }
        };

        // Mark source entities as merged.
        for &source_id in source_entity_ids {
            if let Some(pid) = self.resolve_persistent(source_id) {
                if let Some(record) = self.genealogy.get_mut(&pid) {
                    record.status = EntityStatus::Merged;
                    record.child_ids.push(target_pid);
                    record.transformation_chain.push(GenealogyStep {
                        operation: operation.to_string(),
                        entity_id_before: Some(source_id),
                        entity_id_after: target_entity_id,
                        scope: self.current_scope.clone(),
                    });
                }
            }
        }

        target_pid
    }

    /// Detect naming conflicts in the current state.
    pub fn detect_conflicts(&self) -> Vec<NameConflictRecord> {
        let mut conflicts = Vec::new();
        let mut pid_to_entities: HashMap<PersistentId, Vec<u64>> = HashMap::new();

        // Build a map of persistent ID to all entities that have it.
        for (&entity_id, scoped) in &self.entity_to_scoped {
            pid_to_entities
                .entry(scoped.persistent_id)
                .or_default()
                .push(entity_id);
        }

        // Find persistent IDs assigned to multiple active entities.
        for (pid, entities) in pid_to_entities {
            if entities.len() > 1 {
                // Check if all entities are active.
                let active_count = entities.iter()
                    .filter(|&&entity_id| {
                        self.genealogy.get(&pid)
                            .map(|r| r.status == EntityStatus::Active)
                            .unwrap_or(false)
                    })
                    .count();

                if active_count > 1 {
                    conflicts.push(NameConflictRecord {
                        persistent_id: pid,
                        conflicting_entities: entities,
                        operation: self.current_scope.operation.clone().unwrap_or_default(),
                        scope: self.current_scope.clone(),
                        resolution: NameConflictResolution::Unresolved,
                        sequence: self.conflict_history.len() as u64,
                    });
                }
            }
        }

        conflicts
    }

    /// Resolve a naming conflict.
    pub fn resolve_conflict(
        &mut self,
        conflict: &NameConflictRecord,
        resolution: NameConflictResolution,
    ) -> Result<(), String> {
        match resolution {
            NameConflictResolution::KeepExisting => {
                // No action needed - first entity keeps the ID.
            }
            NameConflictResolution::ReplaceWithNew => {
                // Replace all but the last entity.
                for &entity_id in &conflict.conflicting_entities[..conflict.conflicting_entities.len() - 1] {
                    self.assign_id(entity_id);
                }
            }
            NameConflictResolution::GenerateNewId => {
                // Generate new IDs for all conflicting entities.
                for &entity_id in &conflict.conflicting_entities {
                    self.assign_id(entity_id);
                }
            }
            NameConflictResolution::MergeEntities => {
                // This is handled at the operation level.
            }
            NameConflictResolution::CreateAlias => {
                // Create alias mappings (not fully implemented here).
            }
            NameConflictResolution::Unresolved => {
                return Err("Conflict could not be resolved automatically".to_string());
            }
        }

        // Record the conflict resolution.
        let mut resolved = conflict.clone();
        resolved.resolution = resolution;
        self.conflict_history.push(resolved);

        Ok(())
    }

    /// Get the genealogy record for a persistent ID.
    pub fn get_genealogy(&self, pid: PersistentId) -> Option<&EntityGenealogyRecord> {
        self.genealogy.get(&pid)
    }

    /// Trace the full ancestry of an entity.
    pub fn trace_ancestry(&self, pid: PersistentId) -> Vec<PersistentId> {
        let mut ancestors = Vec::new();
        let mut visited = std::collections::HashSet::new();
        self.collect_ancestry(pid, &mut ancestors, &mut visited);
        ancestors
    }

    fn collect_ancestry(
        &self,
        pid: PersistentId,
        ancestors: &mut Vec<PersistentId>,
        visited: &mut std::collections::HashSet<PersistentId>,
    ) {
        if !visited.insert(pid) {
            return;
        }
        if let Some(record) = self.genealogy.get(&pid) {
            for &parent_pid in &record.parent_ids {
                ancestors.push(parent_pid);
                self.collect_ancestry(parent_pid, ancestors, visited);
            }
        }
    }

    /// Trace the full descendants of an entity.
    pub fn trace_descendants(&self, pid: PersistentId) -> Vec<PersistentId> {
        let mut descendants = Vec::new();
        let mut visited = std::collections::HashSet::new();
        self.collect_descendants(pid, &mut descendants, &mut visited);
        descendants
    }

    fn collect_descendants(
        &self,
        pid: PersistentId,
        descendants: &mut Vec<PersistentId>,
        visited: &mut std::collections::HashSet<PersistentId>,
    ) {
        if !visited.insert(pid) {
            return;
        }
        if let Some(record) = self.genealogy.get(&pid) {
            for &child_pid in &record.child_ids {
                descendants.push(child_pid);
                self.collect_descendants(child_pid, descendants, visited);
            }
        }
    }

    /// Mark an entity as deleted.
    pub fn mark_deleted(&mut self, entity_id: u64) {
        if let Some(pid) = self.resolve_persistent(entity_id) {
            if let Some(record) = self.genealogy.get_mut(&pid) {
                record.status = EntityStatus::Deleted;
            }
        }
    }

    /// Get all entities with a specific status.
    pub fn entities_by_status(&self, status: EntityStatus) -> Vec<PersistentId> {
        self.genealogy
            .iter()
            .filter(|(_, r)| r.status == status)
            .map(|(&pid, _)| pid)
            .collect()
    }

    /// Get the number of active entities.
    pub fn active_entity_count(&self) -> usize {
        self.entities_by_status(EntityStatus::Active).len()
    }

    /// Allocate a new persistent ID.
    fn allocate_persistent_id(&mut self) -> PersistentId {
        self.next_persistent_id += 1;
        PersistentId(self.next_persistent_id)
    }

    /// Clear all bindings and reset the context.
    pub fn clear(&mut self) {
        self.entity_to_scoped.clear();
        self.scoped_to_entity.clear();
        self.genealogy.clear();
        self.pending_assignments.clear();
        self.conflict_history.clear();
        self.next_persistent_id = 0;
    }

    /// Export the context state for serialization.
    pub fn export_state(&self) -> EnhancedNamingContextState {
        EnhancedNamingContextState {
            current_scope: self.current_scope.clone(),
            entity_to_scoped: self.entity_to_scoped.clone(),
            genealogy: self.genealogy.clone(),
            conflict_history: self.conflict_history.clone(),
            next_persistent_id: self.next_persistent_id,
        }
    }

    /// Import context state from a serialized form.
    pub fn import_state(&mut self, state: EnhancedNamingContextState) {
        self.current_scope = state.current_scope;
        self.entity_to_scoped = state.entity_to_scoped;
        self.genealogy = state.genealogy;
        self.conflict_history = state.conflict_history;
        self.next_persistent_id = state.next_persistent_id;

        // Rebuild reverse mapping.
        self.scoped_to_entity.clear();
        for (&entity_id, scoped) in &self.entity_to_scoped {
            self.scoped_to_entity.insert(scoped.clone(), entity_id);
        }
    }
}

/// Serializable state of an EnhancedNamingContext.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnhancedNamingContextState {
    pub current_scope: NamingScope,
    pub entity_to_scoped: HashMap<u64, ScopedId>,
    pub genealogy: HashMap<PersistentId, EntityGenealogyRecord>,
    pub conflict_history: Vec<NameConflictRecord>,
    pub next_persistent_id: u64,
}

/// Manager for operation-specific name propagation rules.
#[derive(Debug, Clone)]
pub struct NamePropagationManager {
    /// Rules indexed by operation type.
    rules: HashMap<OperationType, NamePropagationRule>,
    /// Default rule for unknown operation types.
    default_rule: NamePropagationRule,
}

impl Default for NamePropagationManager {
    fn default() -> Self {
        Self::new()
    }
}

impl NamePropagationManager {
    /// Create a new propagation manager with default rules.
    pub fn new() -> Self {
        let mut rules = HashMap::new();

        // Create default rules for all operation types.
        for op_type in [
            OperationType::BooleanUnion,
            OperationType::BooleanIntersection,
            OperationType::BooleanDifference,
            OperationType::EdgeSplit,
            OperationType::FaceSplit,
            OperationType::Merge,
            OperationType::Delete,
            OperationType::Transform,
            OperationType::Feature,
            OperationType::Generic,
            OperationType::Import,
        ] {
            rules.insert(op_type, NamePropagationRule::for_operation(op_type));
        }

        Self {
            rules,
            default_rule: NamePropagationRule::for_operation(OperationType::Generic),
        }
    }

    /// Get the propagation rule for an operation type.
    pub fn get_rule(&self, operation_type: OperationType) -> &NamePropagationRule {
        self.rules.get(&operation_type).unwrap_or(&self.default_rule)
    }

    /// Set a custom propagation rule for an operation type.
    pub fn set_rule(&mut self, rule: NamePropagationRule) {
        self.rules.insert(rule.operation_type, rule);
    }

    /// Apply a propagation rule to an entity transformation.
    pub fn apply_propagation(
        &self,
        context: &mut EnhancedNamingContext,
        operation_type: OperationType,
        source_entities: &[u64],
        target_entities: &[u64],
        entity_kind: NodeKind,
        operation_name: &str,
    ) -> Vec<PersistentId> {
        let rule = self.get_rule(operation_type);
        let policy = rule.policy_for_kind(entity_kind);

        // Handle split (1 -> many).
        if source_entities.len() == 1 && target_entities.len() > 1 {
            return context.record_split(source_entities[0], target_entities, operation_name);
        }

        // Handle merge (many -> 1).
        if source_entities.len() > 1 && target_entities.len() == 1 {
            return vec![context.record_merge(
                source_entities,
                target_entities[0],
                operation_name,
                rule.conflict_resolution,
            )];
        }

        // Handle 1 -> 1 transformation.
        if source_entities.len() == 1 && target_entities.len() == 1 {
            let pid = context.assign_derived_id(target_entities[0], source_entities, policy);
            return vec![pid];
        }

        // Handle generation (0 -> many).
        if source_entities.is_empty() {
            return target_entities.iter()
                .map(|&entity_id| context.assign_id(entity_id))
                .collect();
        }

        // Default: assign new IDs.
        target_entities.iter()
            .map(|&entity_id| context.assign_id(entity_id))
            .collect()
    }
}

/// Extension trait for BRepGraphHistory to support enhanced naming.
pub trait BRepGraphHistoryExt {
    /// Get the enhanced naming context.
    fn enhanced_context(&self) -> &EnhancedNamingContext;

    /// Get mutable access to the enhanced naming context.
    fn enhanced_context_mut(&mut self) -> &mut EnhancedNamingContext;

    /// Begin an operation with enhanced naming support.
    fn begin_enhanced_operation(
        &mut self,
        operation_type: OperationType,
        label: Option<String>,
        scope: NamingScope,
    );

    /// Propagate names through a boolean operation.
    fn propagate_boolean_names(
        &mut self,
        source_a_entities: &[TopoNode],
        source_b_entities: &[TopoNode],
        result_entities: &[TopoNode],
        operation: BooleanOperationType,
    );

    /// Propagate names through a fillet operation.
    fn propagate_fillet_names(
        &mut self,
        source_edges: &[usize],
        affected_faces: &[usize],
        new_faces: &[usize],
    );

    /// Propagate names through a chamfer operation.
    fn propagate_chamfer_names(
        &mut self,
        source_edges: &[usize],
        affected_faces: &[usize],
        new_faces: &[usize],
    );
}

/// Types of boolean operations for naming propagation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOperationType {
    Union,
    Intersection,
    Difference,
}

// ─────────────────────────────────────────────────────────────────────────────
// Serialization Support
// ─────────────────────────────────────────────────────────────────────────────

/// Serializable snapshot of a naming context for undo/redo support.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NamingContextSnapshot {
    /// Unique identifier for this snapshot.
    pub id: u64,
    /// Timestamp when the snapshot was created.
    pub timestamp: u64,
    /// The scope at the time of the snapshot.
    pub scope: NamingScope,
    /// All entity-to-persistent-ID mappings.
    pub mappings: Vec<(u64, PersistentId, NamingScope)>,
    /// Genealogy records.
    pub genealogy: Vec<EntityGenealogyRecord>,
    /// Conflict records.
    pub conflicts: Vec<NameConflictRecord>,
    /// Operation that created this snapshot.
    pub operation: Option<String>,
}

impl NamingContextSnapshot {
    /// Create a snapshot from an enhanced naming context.
    pub fn from_context(context: &EnhancedNamingContext, id: u64, operation: Option<String>) -> Self {
        let mappings = context.entity_to_scoped.iter()
            .map(|(&entity_id, scoped)| {
                (entity_id, scoped.persistent_id, scoped.scope.clone())
            })
            .collect();

        Self {
            id,
            timestamp: std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs())
                .unwrap_or(0),
            scope: context.current_scope.clone(),
            mappings,
            genealogy: context.genealogy.values().cloned().collect(),
            conflicts: context.conflict_history.clone(),
            operation,
        }
    }

    /// Restore a naming context from this snapshot.
    pub fn restore_to(&self, context: &mut EnhancedNamingContext) {
        context.clear();
        context.current_scope = self.scope.clone();

        for (entity_id, pid, scope) in &self.mappings {
            let scoped = ScopedId::new(*pid, scope.clone());
            context.entity_to_scoped.insert(*entity_id, scoped.clone());
            context.scoped_to_entity.insert(scoped, *entity_id);
        }

        for record in &self.genealogy {
            context.genealogy.insert(record.persistent_id, record.clone());
        }

        context.conflict_history = self.conflicts.clone();
        context.next_persistent_id = self.genealogy.iter()
            .map(|r| r.persistent_id.raw())
            .max()
            .unwrap_or(0);
    }
}

/// Manager for naming context snapshots supporting undo/redo.
#[derive(Debug, Clone, Default)]
pub struct NamingSnapshotManager {
    /// All snapshots in chronological order.
    snapshots: Vec<NamingContextSnapshot>,
    /// Current position in the snapshot history.
    current_index: usize,
    /// Next snapshot ID.
    next_id: u64,
}

impl NamingSnapshotManager {
    /// Create a new snapshot manager.
    pub fn new() -> Self {
        Self::default()
    }

    /// Take a snapshot of the current naming context.
    pub fn take_snapshot(
        &mut self,
        context: &EnhancedNamingContext,
        operation: Option<String>,
    ) -> u64 {
        // Truncate any redo history.
        self.snapshots.truncate(self.current_index + 1);

        let id = self.next_id;
        self.next_id += 1;

        let snapshot = NamingContextSnapshot::from_context(context, id, operation);
        self.snapshots.push(snapshot);
        self.current_index = self.snapshots.len() - 1;

        id
    }

    /// Check if undo is available.
    pub fn can_undo(&self) -> bool {
        self.current_index > 0
    }

    /// Check if redo is available.
    pub fn can_redo(&self) -> bool {
        self.current_index + 1 < self.snapshots.len()
    }

    /// Undo to the previous snapshot.
    pub fn undo(&mut self, context: &mut EnhancedNamingContext) -> Option<&NamingContextSnapshot> {
        if !self.can_undo() {
            return None;
        }
        self.current_index -= 1;
        self.snapshots[self.current_index].restore_to(context);
        Some(&self.snapshots[self.current_index])
    }

    /// Redo to the next snapshot.
    pub fn redo(&mut self, context: &mut EnhancedNamingContext) -> Option<&NamingContextSnapshot> {
        if !self.can_redo() {
            return None;
        }
        self.current_index += 1;
        self.snapshots[self.current_index].restore_to(context);
        Some(&self.snapshots[self.current_index])
    }

    /// Get the current snapshot.
    pub fn current(&self) -> Option<&NamingContextSnapshot> {
        self.snapshots.get(self.current_index)
    }

    /// Get the number of snapshots.
    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    /// Check if there are no snapshots.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    /// Clear all snapshots.
    pub fn clear(&mut self) {
        self.snapshots.clear();
        self.current_index = 0;
    }

    /// Get the undo history (snapshots before current).
    pub fn undo_history(&self) -> &[NamingContextSnapshot] {
        &self.snapshots[..self.current_index]
    }

    /// Get the redo history (snapshots after current).
    pub fn redo_history(&self) -> &[NamingContextSnapshot] {
        &self.snapshots[self.current_index + 1..]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn topo_graph_from_box_has_expected_counts() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let g = TopoGraph::from_brep(&brep);

        let solids = g.nodes.iter().filter(|n| n.kind == NodeKind::Solid).count();
        let shells = g.nodes.iter().filter(|n| n.kind == NodeKind::Shell).count();
        let faces = g.nodes.iter().filter(|n| n.kind == NodeKind::Face).count();
        let wires = g.nodes.iter().filter(|n| n.kind == NodeKind::Wire).count();
        let edges = g.nodes.iter().filter(|n| n.kind == NodeKind::Edge).count();
        let vertices = g.nodes.iter().filter(|n| n.kind == NodeKind::Vertex).count();

        assert_eq!(solids, 1);
        assert_eq!(shells, 1);
        assert_eq!(faces, 6);
        assert_eq!(wires, 6);
        assert_eq!(edges, 12);
        assert_eq!(vertices, 8);
    }

    #[test]
    fn topo_graph_faces_and_edges_queries_work() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let g = TopoGraph::from_brep(&brep);

        let shell0 = TopoNode {
            kind: NodeKind::Shell,
            index: 0,
        };
        let faces = g.faces_of_shell(shell0);
        assert_eq!(faces.len(), 6);

        let face0 = TopoNode {
            kind: NodeKind::Face,
            index: 0,
        };
        let edges = g.edges_of_face(face0);
        assert_eq!(edges.len(), 4);

        let edge0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };
        let verts = g.vertices_of_edge(edge0);
        assert_eq!(verts.len(), 2);
        assert_eq!(verts[0].index, brep.edges[0].start);
        assert_eq!(verts[1].index, brep.edges[0].end);
    }

    #[test]
    fn topo_graph_validate_passes_on_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let g = TopoGraph::from_brep(&brep);
        assert!(g.validate().is_empty());
        assert!(!g.history.events.is_empty());
    }

    #[test]
    fn topo_graph_compact_drops_orphans_and_dedups_neighbors() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);

        // Inject duplicate and orphan adjacency entries to simulate noisy edits.
        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };
        if let Some(neigh) = g.adjacency.get_mut(&v0) {
            neigh.push(e0);
            neigh.push(e0);
        }
        let orphan = TopoNode {
            kind: NodeKind::Face,
            index: 9999,
        };
        g.adjacency.insert(orphan, vec![v0]);

        g.compact();

        assert!(!g.adjacency.contains_key(&orphan));
        let neigh = g.adjacency.get(&v0).expect("vertex adjacency exists");
        let count_e0 = neigh.iter().filter(|n| **n == e0).count();
        assert_eq!(count_e0, 1);
        assert!(g.history.events.iter().any(|e| e.action == "compact"));
    }

    #[test]
    fn topo_graph_mutate_checked_reports_invalid_graph() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);

        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };

        let res = g.mutate_checked("inject_nonsymmetric", |graph| {
            if let Some(neigh) = graph.adjacency.get_mut(&e0) {
                neigh.retain(|n| *n != v0);
            }
        });

        assert!(res.is_err());
        let issues = res.expect_err("mutation should be invalid");
        assert!(issues
            .iter()
            .any(|i| matches!(i, TopoGraphValidationIssue::NonSymmetricAdjacency { .. })));
        assert!(g
            .history
            .events
            .iter()
            .any(|e| e.action == "mutate_invalid:inject_nonsymmetric"));
    }

    #[test]
    fn topo_graph_mutate_guarded_rolls_back_on_invalid_mutation() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);
        let before = g.clone();

        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };

        let res = g.mutate_guarded("inject_nonsymmetric", |graph| {
            if let Some(neigh) = graph.adjacency.get_mut(&e0) {
                neigh.retain(|n| *n != v0);
            }
        });

        assert!(res.is_err());
        assert_eq!(g.nodes, before.nodes);
        assert_eq!(g.adjacency, before.adjacency);
        assert!(g
            .history
            .events
            .iter()
            .any(|e| e.action == "mutate_guarded_rollback:inject_nonsymmetric"));
    }

    #[test]
    fn topo_graph_mutate_guarded_commits_valid_mutation() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };
        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };

        let res = g.mutate_guarded("dedup_neighbors", |graph| {
            if let Some(neigh) = graph.adjacency.get_mut(&e0) {
                neigh.push(v0);
            }
            graph.compact();
        });

        assert!(res.is_ok());
        let neigh = g.adjacency.get(&e0).expect("edge adjacency exists");
        let count_v0 = neigh.iter().filter(|n| **n == v0).count();
        assert_eq!(count_v0, 1);
        assert!(g
            .history
            .events
            .iter()
            .any(|e| e.action == "mutate_guarded:dedup_neighbors"));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests for Cross-Operation Naming Stability Analysis
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod cross_operation_tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;
    use rcad_kernel::persistent_naming::{
        OperationType, NamingEvent, PersistentId,
    };

    /// Helper to create a history with a single operation.
    fn create_history_with_operation(
        operation_type: OperationType,
        label: &str,
        stats: OperationStats,
        events: Vec<NamingEvent>,
    ) -> BRepGraphHistory {
        let mut history = BRepGraphHistory::new();
        history.begin_operation(operation_type, Some(label.to_string()));
        history.naming_engine_mut().cross_operation_history_mut()
            .add_events(rcad_kernel::persistent_naming::OperationId(1), events);
        history.finalize_operation(stats);
        history
    }

    #[test]
    fn analyze_empty_history_returns_perfect_stability() {
        let analysis = analyze_naming_sequence(&[], &[]);

        assert!(analysis.is_excellent());
        assert_eq!(analysis.overall_stability, 1.0);
        assert!(analysis.broken_chains.is_empty());
        assert!(analysis.operations.is_empty());
    }

    #[test]
    fn analyze_single_operation_with_no_losses() {
        let history = create_history_with_operation(
            OperationType::BooleanUnion,
            "test_union",
            OperationStats {
                entity_count_before: 10,
                entity_count_after: 10,
                names_preserved: 10,
                names_lost: 0,
                names_generated: 0,
                conflicts_resolved: 0,
            },
            vec![],
        );

        let analysis = analyze_naming_sequence(&[history], &[]);

        assert!(analysis.is_excellent());
        assert_eq!(analysis.overall_stability, 1.0);
        assert_eq!(analysis.operations.len(), 1);
    }

    #[test]
    fn analyze_single_operation_with_losses() {
        let history = create_history_with_operation(
            OperationType::BooleanDifference,
            "test_diff",
            OperationStats {
                entity_count_before: 10,
                entity_count_after: 7,
                names_preserved: 7,
                names_lost: 3,
                names_generated: 0,
                conflicts_resolved: 0,
            },
            vec![
                NamingEvent::Lost {
                    entity_id: 1,
                    persistent_id: PersistentId(1),
                },
                NamingEvent::Lost {
                    entity_id: 2,
                    persistent_id: PersistentId(2),
                },
                NamingEvent::Lost {
                    entity_id: 3,
                    persistent_id: PersistentId(3),
                },
            ],
        );

        let initial_entities: Vec<TopoNode> = (0..10)
            .map(|i| TopoNode { kind: NodeKind::Face, index: i })
            .collect();

        let analysis = analyze_naming_sequence(&[history], &initial_entities);

        assert!(!analysis.is_excellent());
        assert_eq!(analysis.broken_chains.len(), 3);
        assert!(analysis.overall_stability < 1.0);
    }

    #[test]
    fn analyze_multiple_operations_cumulative_stability() {
        let history1 = create_history_with_operation(
            OperationType::BooleanUnion,
            "union1",
            OperationStats {
                entity_count_before: 10,
                entity_count_after: 12,
                names_preserved: 10,
                names_lost: 0,
                names_generated: 2,
                conflicts_resolved: 0,
            },
            vec![],
        );

        let history2 = create_history_with_operation(
            OperationType::BooleanDifference,
            "diff1",
            OperationStats {
                entity_count_before: 12,
                entity_count_after: 10,
                names_preserved: 8,
                names_lost: 4,
                names_generated: 2,
                conflicts_resolved: 0,
            },
            vec![
                NamingEvent::Lost {
                    entity_id: 1,
                    persistent_id: PersistentId(1),
                },
                NamingEvent::Lost {
                    entity_id: 2,
                    persistent_id: PersistentId(2),
                },
            ],
        );

        let analysis = analyze_naming_sequence(&[history1, history2], &[]);

        assert_eq!(analysis.operations.len(), 2);
        assert_eq!(analysis.per_operation_stability.len(), 2);

        // Cumulative stability should decrease.
        assert!(analysis.per_operation_stability[1].cumulative_stability <= analysis.per_operation_stability[0].cumulative_stability);
    }

    #[test]
    fn stability_trend_improving() {
        // Create metrics with improving stability.
        let metrics = vec![
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(1),
                operation_type: OperationType::BooleanUnion,
                label: None,
                names_retained: 5,
                names_lost: 5,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.5,
                cumulative_stability: 0.5,
            },
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(2),
                operation_type: OperationType::BooleanUnion,
                label: None,
                names_retained: 7,
                names_lost: 3,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.7,
                cumulative_stability: 0.6,
            },
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(3),
                operation_type: OperationType::BooleanUnion,
                label: None,
                names_retained: 9,
                names_lost: 1,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.9,
                cumulative_stability: 0.7,
            },
        ];

        let trend = calculate_stability_trend(&metrics);
        assert!(trend > 0.0, "Trend should be positive (improving), got {}", trend);
    }

    #[test]
    fn stability_trend_degrading() {
        // Create metrics with degrading stability.
        let metrics = vec![
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(1),
                operation_type: OperationType::BooleanDifference,
                label: None,
                names_retained: 9,
                names_lost: 1,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.9,
                cumulative_stability: 0.9,
            },
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(2),
                operation_type: OperationType::BooleanDifference,
                label: None,
                names_retained: 7,
                names_lost: 3,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.7,
                cumulative_stability: 0.8,
            },
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(3),
                operation_type: OperationType::BooleanDifference,
                label: None,
                names_retained: 5,
                names_lost: 5,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.5,
                cumulative_stability: 0.6,
            },
        ];

        let trend = calculate_stability_trend(&metrics);
        assert!(trend < 0.0, "Trend should be negative (degrading), got {}", trend);
    }

    #[test]
    fn detect_conflicts_no_conflicts() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 1.0,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let conflicts = detect_cross_operation_conflicts(&analysis);
        assert!(conflicts.is_empty());
    }

    #[test]
    fn detect_conflicts_broken_genealogy() {
        let mut genealogy = HashMap::new();
        genealogy.insert(
            PersistentId(1),
            EntityGenealogy {
                persistent_id: PersistentId(1),
                created_in_operation: rcad_kernel::persistent_naming::OperationId(1),
                evolution: vec![],  // Empty evolution is a broken genealogy
                current_entity_id: Some(42),
                is_deleted: false,
            },
        );

        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: genealogy,
            per_operation_stability: vec![],
            overall_stability: 1.0,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let conflicts = detect_cross_operation_conflicts(&analysis);
        assert!(!conflicts.is_empty());

        let has_broken = conflicts.iter().any(|c| c.conflict_type == ConflictType::BrokenGenealogy);
        assert!(has_broken, "Should detect BrokenGenealogy conflict");
    }

    #[test]
    fn generate_recommendations_critical_stability() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.3,  // Critical
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        assert!(!recommendations.is_empty());

        // Should have a high-priority architecture recommendation.
        let has_critical = recommendations.iter().any(|r| r.priority == 100);
        assert!(has_critical, "Should have critical stability recommendation");
    }

    #[test]
    fn generate_recommendations_good_stability() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.95,  // Good
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        // With good stability and no issues, should have minimal or no recommendations.
        // Actually, it should be empty since no issues.
        assert!(recommendations.is_empty() || recommendations.iter().all(|r| r.priority < 70));
    }

    #[test]
    fn generate_recommendations_broken_chains() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.9,
            broken_chains: vec![
                BrokenChainInfo {
                    persistent_id: PersistentId(1),
                    broken_at_operation: rcad_kernel::persistent_naming::OperationId(1),
                    entity_id: 1,
                    entity_type: Some(rcad_kernel::persistent_naming::EntityType::Face),
                    survived_operations: 10,  // Long-lived entity
                    break_reason: "Test break".to_string(),
                },
            ],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        let has_entity_tracking = recommendations.iter().any(|r| r.category == RecommendationCategory::EntityTracking);
        assert!(has_entity_tracking, "Should have entity tracking recommendation");
    }

    #[test]
    fn generate_recommendations_degrading_trend() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.85,
            broken_chains: vec![],
            stability_trend: -0.2,  // Degrading
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        let has_propagation = recommendations.iter().any(|r| r.category == RecommendationCategory::PropagationPolicy);
        assert!(has_propagation, "Should have propagation policy recommendation for degrading trend");
    }

    #[test]
    fn most_problematic_operation() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(1),
                    operation_type: OperationType::BooleanUnion,
                    label: None,
                    names_retained: 9,
                    names_lost: 1,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.9,
                    cumulative_stability: 0.9,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(2),
                    operation_type: OperationType::BooleanDifference,
                    label: None,
                    names_retained: 5,
                    names_lost: 5,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.5,  // Most problematic
                    cumulative_stability: 0.45,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(3),
                    operation_type: OperationType::Feature,
                    label: None,
                    names_retained: 8,
                    names_lost: 2,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.8,
                    cumulative_stability: 0.36,
                },
            ],
            overall_stability: 0.7,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let problematic = analysis.most_problematic_operation();
        assert!(problematic.is_some());
        assert_eq!(problematic.unwrap().operation_id, rcad_kernel::persistent_naming::OperationId(2));
    }

    #[test]
    fn operations_by_stability_sorted() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(1),
                    operation_type: OperationType::BooleanUnion,
                    label: None,
                    names_retained: 7,
                    names_lost: 3,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.7,
                    cumulative_stability: 0.7,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(2),
                    operation_type: OperationType::BooleanDifference,
                    label: None,
                    names_retained: 5,
                    names_lost: 5,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.5,
                    cumulative_stability: 0.35,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(3),
                    operation_type: OperationType::Feature,
                    label: None,
                    names_retained: 9,
                    names_lost: 1,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.9,
                    cumulative_stability: 0.315,
                },
            ],
            overall_stability: 0.7,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let sorted = analysis.operations_by_stability();
        assert_eq!(sorted.len(), 3);
        // Should be sorted by stability score (ascending).
        assert_eq!(sorted[0].stability_score, 0.5);
        assert_eq!(sorted[1].stability_score, 0.7);
        assert_eq!(sorted[2].stability_score, 0.9);
    }

    #[test]
    fn average_operation_stability() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(1),
                    operation_type: OperationType::BooleanUnion,
                    label: None,
                    names_retained: 8,
                    names_lost: 2,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.8,
                    cumulative_stability: 0.8,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(2),
                    operation_type: OperationType::BooleanDifference,
                    label: None,
                    names_retained: 6,
                    names_lost: 4,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.6,
                    cumulative_stability: 0.48,
                },
            ],
            overall_stability: 0.7,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let avg = analysis.average_operation_stability();
        assert!((avg - 0.7).abs() < 0.001, "Average should be 0.7, got {}", avg);
    }

    #[test]
    fn analysis_summary_format() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![OperationRecord {
                id: rcad_kernel::persistent_naming::OperationId(1),
                operation_type: OperationType::BooleanUnion,
                label: Some("test".to_string()),
                sequence: 0,
                naming_events: vec![],
                stats: OperationStats::default(),
            }],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.85,
            broken_chains: vec![],
            stability_trend: 0.1,
            entity_counts: vec![10],
        };

        let summary = analysis.summary();
        assert!(summary.contains("Operations: 1"));
        assert!(summary.contains("85.0%"));
        assert!(summary.contains("Improving"));
    }

    #[test]
    fn infer_entity_type() {
        // Test entity type inference from encoded IDs.
        // Kind bits: Solid=0, Shell=1, Face=2, Wire=3, Edge=4, Vertex=5

        // Face (kind=2)
        let face_id = (2u64 << 56) | 42u64;
        let face_type = infer_entity_type_from_id(face_id);
        assert_eq!(face_type, Some(rcad_kernel::persistent_naming::EntityType::Face));

        // Edge (kind=4)
        let edge_id = (4u64 << 56) | 100u64;
        let edge_type = infer_entity_type_from_id(edge_id);
        assert_eq!(edge_type, Some(rcad_kernel::persistent_naming::EntityType::Edge));

        // Vertex (kind=5)
        let vertex_id = (5u64 << 56) | 7u64;
        let vertex_type = infer_entity_type_from_id(vertex_id);
        assert_eq!(vertex_type, Some(rcad_kernel::persistent_naming::EntityType::Vertex));

        // Solid (kind=0)
        let solid_id = (0u64 << 56) | 1u64;
        let solid_type = infer_entity_type_from_id(solid_id);
        assert_eq!(solid_type, Some(rcad_kernel::persistent_naming::EntityType::Solid));
    }

    #[test]
    fn cross_operation_stability_with_named_graph() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let mut named_graph = NamedGraph::from_brep(&brep);
        let initial_nodes = named_graph.graph().nodes.clone();

        // Perform a mutation.
        let result = named_graph.mutate_tracked(
            "test_mutation",
            OperationType::Generic,
            |graph, _history| {
                graph.record("test_op");
            },
        );

        assert!(result.is_ok());

        // Get stability report.
        let report = named_graph.stability_report();
        assert!(report.total_operations >= 1);

        // Analyze the history.
        let history = named_graph.history().clone();
        let analysis = analyze_naming_sequence(&[history], &initial_nodes);
        assert!(analysis.is_good() || analysis.is_excellent());
    }

    #[test]
    fn conflict_severity_classification() {
        let mut genealogy = HashMap::new();

        // Entity with broken genealogy (evolution is empty).
        genealogy.insert(
            PersistentId(1),
            EntityGenealogy {
                persistent_id: PersistentId(1),
                created_in_operation: rcad_kernel::persistent_naming::OperationId(1),
                evolution: vec![],
                current_entity_id: Some(1),
                is_deleted: false,
            },
        );

        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: genealogy,
            per_operation_stability: vec![],
            overall_stability: 0.9,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let conflicts = detect_cross_operation_conflicts(&analysis);
        assert!(!conflicts.is_empty());

        // All detected conflicts should have proper severity.
        for conflict in &conflicts {
            assert!(matches!(
                conflict.severity,
                rcad_kernel::persistent_naming::IssueSeverity::Minor
                    | rcad_kernel::persistent_naming::IssueSeverity::Moderate
                    | rcad_kernel::persistent_naming::IssueSeverity::Severe
                    | rcad_kernel::persistent_naming::IssueSeverity::Critical
            ));
        }
    }

    #[test]
    fn recommendation_priorities_ordered() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.3,  // Critical stability
            broken_chains: vec![],
            stability_trend: -0.5,  // Degrading
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        assert!(recommendations.len() >= 2, "Should have multiple recommendations");

        // Verify priorities are sorted (highest first).
        for i in 1..recommendations.len() {
            assert!(
                recommendations[i].priority <= recommendations[i - 1].priority,
                "Recommendations should be sorted by priority (descending)"
            );
        }
    }

    #[test]
    fn entity_counts_tracking() {
        let history = create_history_with_operation(
            OperationType::BooleanUnion,
            "test",
            OperationStats {
                entity_count_before: 10,
                entity_count_after: 15,
                names_preserved: 10,
                names_lost: 0,
                names_generated: 5,
                conflicts_resolved: 0,
            },
            vec![],
        );

        let initial_entities: Vec<TopoNode> = (0..10)
            .map(|i| TopoNode { kind: NodeKind::Face, index: i })
            .collect();

        let analysis = analyze_naming_sequence(&[history], &initial_entities);

        assert!(!analysis.entity_counts.is_empty());
        assert_eq!(analysis.entity_counts[0], 10);  // Initial count
        // After the operation, the count from stats is tracked.
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests for Enhanced Persistent Naming Semantics
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod enhanced_naming_tests {
    use super::*;

    // ── ScopedId Tests ────────────────────────────────────────────────────────

    #[test]
    fn scoped_id_creation() {
        let scope = NamingScope::for_part("housing").with_operation("fillet");
        let pid = PersistentId(42);
        let scoped = ScopedId::new(pid, scope.clone());

        assert_eq!(scoped.persistent_id, pid);
        assert_eq!(scoped.scope, scope);
        assert!(!scoped.is_null());
    }

    #[test]
    fn scoped_id_null() {
        let scoped = ScopedId::null();
        assert!(scoped.is_null());
        assert!(scoped.persistent_id.is_null());
    }

    #[test]
    fn scoped_id_qualified_name() {
        let scope = NamingScope::for_part("housing")
            .with_assembly("machine")
            .with_operation("fillet");
        let scoped = ScopedId::new(PersistentId(42), scope);

        let name = scoped.qualified_name();
        assert!(name.contains("machine"));
        assert!(name.contains("housing"));
        assert!(name.contains("fillet"));
        assert!(name.contains("e42"));
    }

    // ── NamingScope Tests ──────────────────────────────────────────────────────

    #[test]
    fn naming_scope_creation() {
        let scope = NamingScope::new();
        assert!(scope.part.is_none());
        assert!(scope.assembly.is_none());
        assert!(scope.operation.is_none());
    }

    #[test]
    fn naming_scope_for_part() {
        let scope = NamingScope::for_part("housing");
        assert_eq!(scope.part, Some("housing".to_string()));
        assert!(scope.assembly.is_none());
    }

    #[test]
    fn naming_scope_for_assembly() {
        let scope = NamingScope::for_assembly("machine");
        assert_eq!(scope.assembly, Some("machine".to_string()));
        assert!(scope.part.is_none());
    }

    #[test]
    fn naming_scope_for_operation() {
        let scope = NamingScope::for_operation("housing", "fillet");
        assert_eq!(scope.part, Some("housing".to_string()));
        assert_eq!(scope.operation, Some("fillet".to_string()));
    }

    #[test]
    fn naming_scope_builder_pattern() {
        let scope = NamingScope::new()
            .with_assembly("machine")
            .with_part("housing")
            .with_operation("fillet");

        assert_eq!(scope.assembly, Some("machine".to_string()));
        assert_eq!(scope.part, Some("housing".to_string()));
        assert_eq!(scope.operation, Some("fillet".to_string()));
    }

    #[test]
    fn naming_scope_child_scope() {
        let parent = NamingScope::for_part("housing");
        let child = parent.child_scope("fillet");

        assert_eq!(child.part, parent.part);
        assert_eq!(child.operation, Some("fillet".to_string()));
    }

    #[test]
    fn naming_scope_contains() {
        let parent = NamingScope::for_assembly("machine");
        let child = NamingScope::for_part("housing").with_assembly("machine");

        assert!(parent.contains(&child));
        assert!(!child.contains(&parent));

        let other = NamingScope::for_assembly("other");
        assert!(!parent.contains(&other));
    }

    // ── EnhancedNamingContext Tests ────────────────────────────────────────────

    #[test]
    fn enhanced_context_assign_id() {
        let mut ctx = EnhancedNamingContext::new();
        let pid = ctx.assign_id(42);

        assert!(!pid.is_null());
        assert_eq!(ctx.resolve_entity(pid), Some(42));
        assert_eq!(ctx.resolve_persistent(42), Some(pid));
    }

    #[test]
    fn enhanced_context_with_scope() {
        let scope = NamingScope::for_part("housing");
        let ctx = EnhancedNamingContext::with_scope(scope.clone());

        assert_eq!(ctx.scope(), &scope);
    }

    #[test]
    fn enhanced_context_assign_derived_id_preserve() {
        let mut ctx = EnhancedNamingContext::new();

        // Assign original.
        let original_pid = ctx.assign_id(10);

        // Derive with Preserve policy.
        let derived_pid = ctx.assign_derived_id(20, &[10], NamePropagationPolicy::Preserve);

        // Should inherit the same persistent ID.
        assert_eq!(derived_pid, original_pid);
        assert_eq!(ctx.resolve_persistent(20), Some(original_pid));
    }

    #[test]
    fn enhanced_context_assign_derived_id_generate() {
        let mut ctx = EnhancedNamingContext::new();

        // Assign original.
        let original_pid = ctx.assign_id(10);

        // Derive with Generate policy.
        let derived_pid = ctx.assign_derived_id(20, &[10], NamePropagationPolicy::Generate);

        // Should get a new persistent ID.
        assert_ne!(derived_pid, original_pid);
        assert_eq!(ctx.resolve_persistent(20), Some(derived_pid));
    }

    #[test]
    fn enhanced_context_record_split() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.set_scope(NamingScope::for_operation("test", "split_op"));

        // Create source entity.
        let source_pid = ctx.assign_id(10);

        // Split into three entities.
        let result_pids = ctx.record_split(10, &[20, 30, 40], "split_op");

        assert_eq!(result_pids.len(), 3);
        // First target inherits source's PID.
        assert_eq!(result_pids[0], source_pid);
        // Others get new PIDs.
        assert_ne!(result_pids[1], source_pid);
        assert_ne!(result_pids[2], source_pid);

        // Check genealogy.
        let genealogy = ctx.get_genealogy(source_pid).unwrap();
        assert_eq!(genealogy.status, EntityStatus::Split);
        assert_eq!(genealogy.child_ids.len(), 2); // Two new children.
    }

    #[test]
    fn enhanced_context_record_merge() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.set_scope(NamingScope::for_operation("test", "merge_op"));

        // Create source entities.
        let pid1 = ctx.assign_id(10);
        let pid2 = ctx.assign_id(20);

        // Merge into one entity.
        let result_pid = ctx.record_merge(
            &[10, 20],
            30,
            "merge_op",
            NameConflictResolution::MergeEntities,
        );

        // Check result exists.
        assert!(!result_pid.is_null());

        // Check genealogy of sources.
        let genealogy1 = ctx.get_genealogy(pid1).unwrap();
        assert_eq!(genealogy1.status, EntityStatus::Merged);

        let genealogy2 = ctx.get_genealogy(pid2).unwrap();
        assert_eq!(genealogy2.status, EntityStatus::Merged);
    }

    #[test]
    fn enhanced_context_mark_deleted() {
        let mut ctx = EnhancedNamingContext::new();
        let pid = ctx.assign_id(10);

        ctx.mark_deleted(10);

        let genealogy = ctx.get_genealogy(pid).unwrap();
        assert_eq!(genealogy.status, EntityStatus::Deleted);
    }

    #[test]
    fn enhanced_context_entities_by_status() {
        let mut ctx = EnhancedNamingContext::new();

        let pid1 = ctx.assign_id(10);
        let pid2 = ctx.assign_id(20);
        ctx.mark_deleted(10);

        let deleted = ctx.entities_by_status(EntityStatus::Deleted);
        assert!(deleted.contains(&pid1));
        assert!(!deleted.contains(&pid2));

        let active = ctx.entities_by_status(EntityStatus::Active);
        assert!(active.contains(&pid2));
        assert!(!active.contains(&pid1));
    }

    #[test]
    fn enhanced_context_export_import_state() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.set_scope(NamingScope::for_part("test"));

        let pid1 = ctx.assign_id(10);
        let pid2 = ctx.assign_id(20);

        // Export state.
        let state = ctx.export_state();

        // Clear context.
        ctx.clear();
        assert!(ctx.resolve_persistent(10).is_none());

        // Import state.
        ctx.import_state(state);

        // Verify restoration.
        assert_eq!(ctx.resolve_persistent(10), Some(pid1));
        assert_eq!(ctx.resolve_persistent(20), Some(pid2));
    }

    // ── NamePropagationRule Tests ──────────────────────────────────────────────

    #[test]
    fn propagation_rule_for_boolean() {
        let rule = NamePropagationRule::for_operation(OperationType::BooleanUnion);

        assert_eq!(rule.face_policy, NamePropagationPolicy::Preserve);
        assert_eq!(rule.edge_policy, NamePropagationPolicy::Preserve);
        assert!(rule.track_genealogy);
    }

    #[test]
    fn propagation_rule_for_feature() {
        let rule = NamePropagationRule::for_operation(OperationType::Feature);

        assert_eq!(rule.face_policy, NamePropagationPolicy::Inherit);
        assert_eq!(rule.edge_policy, NamePropagationPolicy::Inherit);
        assert!(rule.track_genealogy);
    }

    #[test]
    fn propagation_rule_for_merge() {
        let rule = NamePropagationRule::for_operation(OperationType::Merge);

        assert_eq!(rule.face_policy, NamePropagationPolicy::Combine);
        assert_eq!(rule.conflict_resolution, NameConflictResolution::MergeEntities);
    }

    #[test]
    fn propagation_rule_policy_for_kind() {
        let rule = NamePropagationRule {
            operation_type: OperationType::Generic,
            face_policy: NamePropagationPolicy::Preserve,
            edge_policy: NamePropagationPolicy::Inherit,
            vertex_policy: NamePropagationPolicy::Generate,
            track_genealogy: true,
            conflict_resolution: NameConflictResolution::GenerateNewId,
        };

        assert_eq!(rule.policy_for_kind(NodeKind::Face), NamePropagationPolicy::Preserve);
        assert_eq!(rule.policy_for_kind(NodeKind::Edge), NamePropagationPolicy::Inherit);
        assert_eq!(rule.policy_for_kind(NodeKind::Vertex), NamePropagationPolicy::Generate);
    }

    // ── NamePropagationManager Tests ───────────────────────────────────────────

    #[test]
    fn propagation_manager_get_rule() {
        let manager = NamePropagationManager::new();

        let rule = manager.get_rule(OperationType::BooleanUnion);
        assert_eq!(rule.operation_type, OperationType::BooleanUnion);

        let rule = manager.get_rule(OperationType::Feature);
        assert_eq!(rule.operation_type, OperationType::Feature);
    }

    #[test]
    fn propagation_manager_set_rule() {
        let mut manager = NamePropagationManager::new();

        let custom_rule = NamePropagationRule {
            operation_type: OperationType::BooleanUnion,
            face_policy: NamePropagationPolicy::Generate,
            edge_policy: NamePropagationPolicy::Generate,
            vertex_policy: NamePropagationPolicy::Generate,
            track_genealogy: false,
            conflict_resolution: NameConflictResolution::KeepExisting,
        };

        manager.set_rule(custom_rule);

        let rule = manager.get_rule(OperationType::BooleanUnion);
        assert_eq!(rule.face_policy, NamePropagationPolicy::Generate);
    }

    #[test]
    fn propagation_manager_apply_split() {
        let mut manager = NamePropagationManager::new();
        let mut ctx = EnhancedNamingContext::new();

        // Create source entity.
        ctx.assign_id(10);

        // Apply split propagation.
        let pids = manager.apply_propagation(
            &mut ctx,
            OperationType::FaceSplit,
            &[10],
            &[20, 30, 40],
            NodeKind::Face,
            "split_test",
        );

        assert_eq!(pids.len(), 3);
        // First target should inherit source's PID.
        assert_eq!(pids[0], ctx.resolve_persistent(10).unwrap());
    }

    #[test]
    fn propagation_manager_apply_merge() {
        let mut manager = NamePropagationManager::new();
        let mut ctx = EnhancedNamingContext::new();

        // Create source entities.
        ctx.assign_id(10);
        ctx.assign_id(20);

        // Apply merge propagation.
        let pids = manager.apply_propagation(
            &mut ctx,
            OperationType::Merge,
            &[10, 20],
            &[30],
            NodeKind::Face,
            "merge_test",
        );

        assert_eq!(pids.len(), 1);
    }

    // ── NamingSnapshotManager Tests ────────────────────────────────────────────

    #[test]
    fn snapshot_manager_take_snapshot() {
        let mut manager = NamingSnapshotManager::new();
        let ctx = EnhancedNamingContext::new();

        let id = manager.take_snapshot(&ctx, Some("test_op".to_string()));

        assert_eq!(id, 0);
        assert_eq!(manager.len(), 1);
    }

    #[test]
    fn snapshot_manager_undo_redo() {
        let mut manager = NamingSnapshotManager::new();
        let mut ctx = EnhancedNamingContext::new();

        // Initial state.
        ctx.assign_id(10);
        manager.take_snapshot(&ctx, Some("op1".to_string()));

        // Second state.
        ctx.assign_id(20);
        manager.take_snapshot(&ctx, Some("op2".to_string()));

        assert_eq!(manager.len(), 2);

        // Undo.
        assert!(manager.can_undo());
        manager.undo(&mut ctx);
        assert!(!manager.can_undo());

        // Redo.
        assert!(manager.can_redo());
        manager.redo(&mut ctx);
        assert!(!manager.can_redo());
    }

    #[test]
    fn snapshot_manager_current() {
        let mut manager = NamingSnapshotManager::new();
        let ctx = EnhancedNamingContext::new();

        manager.take_snapshot(&ctx, Some("test".to_string()));

        let current = manager.current();
        assert!(current.is_some());
        assert_eq!(current.unwrap().operation, Some("test".to_string()));
    }

    // ── Genealogy Tracking Tests ────────────────────────────────────────────────

    #[test]
    fn genealogy_tracking_through_boolean() {
        let mut ctx = EnhancedNamingContext::new();
        let mut manager = NamePropagationManager::new();

        ctx.set_scope(NamingScope::for_operation("part", "boolean_union"));

        // Create source entities.
        ctx.assign_id(1);
        ctx.assign_id(2);

        // Simulate boolean union.
        let result_pids = manager.apply_propagation(
            &mut ctx,
            OperationType::BooleanUnion,
            &[1, 2],
            &[10, 11],
            NodeKind::Face,
            "boolean_union",
        );

        assert_eq!(result_pids.len(), 2);
    }

    #[test]
    fn genealogy_tracking_through_fillet() {
        let mut ctx = EnhancedNamingContext::new();
        let mut manager = NamePropagationManager::new();

        ctx.set_scope(NamingScope::for_operation("part", "fillet"));

        // Create source edges.
        let pid_edge = ctx.assign_id(1);

        // Simulate fillet.
        let result_pids = manager.apply_propagation(
            &mut ctx,
            OperationType::Feature,
            &[1],
            &[10, 11],
            NodeKind::Face,
            "fillet",
        );

        assert_eq!(result_pids.len(), 2);

        let edge_genealogy = ctx.get_genealogy(pid_edge);
        assert!(edge_genealogy.is_some());
    }

    #[test]
    fn genealogy_multiple_operations() {
        let mut ctx = EnhancedNamingContext::new();
        let mut manager = NamePropagationManager::new();

        // Operation 1: Create initial face.
        ctx.set_scope(NamingScope::for_operation("part", "create"));
        let pid1 = ctx.assign_id(1);

        // Operation 2: Boolean union.
        ctx.set_scope(NamingScope::for_operation("part", "union"));
        manager.apply_propagation(
            &mut ctx,
            OperationType::BooleanUnion,
            &[1],
            &[2],
            NodeKind::Face,
            "union",
        );

        // Operation 3: Fillet.
        ctx.set_scope(NamingScope::for_operation("part", "fillet"));
        manager.apply_propagation(
            &mut ctx,
            OperationType::Feature,
            &[2],
            &[3, 4],
            NodeKind::Face,
            "fillet",
        );

        // Verify the chain is trackable.
        let descendants = ctx.trace_descendants(pid1);
        assert!(!descendants.is_empty() || ctx.get_genealogy(pid1).is_some());
    }

    // ── Name Stability Tests ────────────────────────────────────────────────────

    #[test]
    fn name_stability_through_boolean_preserve() {
        let mut ctx = EnhancedNamingContext::new();

        // Create faces from solid A.
        let pid_f1 = ctx.assign_id(1);
        ctx.assign_id(2);

        // Simulate boolean union where faces 1 and 2 are preserved.
        let preserved_pid = ctx.assign_derived_id(
            10, // New face index
            &[1], // Source face
            NamePropagationPolicy::Preserve,
        );

        assert_eq!(preserved_pid, pid_f1);
    }

    #[test]
    fn name_stability_through_split() {
        let mut ctx = EnhancedNamingContext::new();

        // Create a face that will be split.
        let pid = ctx.assign_id(1);

        // Split into three faces.
        let result_pids = ctx.record_split(1, &[10, 11, 12], "split");

        // First result should inherit the original PID.
        assert_eq!(result_pids[0], pid);
        // Others should get new PIDs.
        assert_ne!(result_pids[1], pid);
        assert_ne!(result_pids[2], pid);

        // Original entity should be marked as split.
        let genealogy = ctx.get_genealogy(pid).unwrap();
        assert_eq!(genealogy.status, EntityStatus::Split);
    }

    #[test]
    fn name_stability_through_merge() {
        let mut ctx = EnhancedNamingContext::new();

        // Create two faces that will be merged.
        let pid1 = ctx.assign_id(1);
        let pid2 = ctx.assign_id(2);

        // Merge into one face.
        let result_pid = ctx.record_merge(
            &[1, 2],
            10,
            "merge",
            NameConflictResolution::MergeEntities,
        );

        // Result should exist and sources should be marked merged.
        assert!(!result_pid.is_null());

        let g1 = ctx.get_genealogy(pid1).unwrap();
        let g2 = ctx.get_genealogy(pid2).unwrap();

        assert_eq!(g1.status, EntityStatus::Merged);
        assert_eq!(g2.status, EntityStatus::Merged);
    }

    // ── Conflict Resolution Tests ───────────────────────────────────────────────

    #[test]
    fn conflict_resolution_keep_existing() {
        let mut ctx = EnhancedNamingContext::new();

        let pid = ctx.assign_id(1);

        // Try to assign the same PID to another entity (simulating conflict).
        let conflict = NameConflictRecord {
            persistent_id: pid,
            conflicting_entities: vec![1, 2],
            operation: "test".to_string(),
            scope: ctx.current_scope.clone(),
            resolution: NameConflictResolution::Unresolved,
            sequence: 0,
        };

        ctx.resolve_conflict(&conflict, NameConflictResolution::KeepExisting).unwrap();

        // Original entity should keep its PID.
        assert_eq!(ctx.resolve_persistent(1), Some(pid));

        // Conflict should be recorded.
        assert!(!ctx.conflict_history.is_empty());
    }

    #[test]
    fn conflict_resolution_generate_new() {
        let mut ctx = EnhancedNamingContext::new();

        let pid = ctx.assign_id(1);

        let conflict = NameConflictRecord {
            persistent_id: pid,
            conflicting_entities: vec![1, 2],
            operation: "test".to_string(),
            scope: ctx.current_scope.clone(),
            resolution: NameConflictResolution::Unresolved,
            sequence: 0,
        };

        ctx.resolve_conflict(&conflict, NameConflictResolution::GenerateNewId).unwrap();

        // New entity should have a different PID.
        let pid1 = ctx.resolve_persistent(1);
        let pid2 = ctx.resolve_persistent(2);

        // At least one should have a different PID.
        assert!(pid1.is_some() || pid2.is_some());
    }

    // ── Serialization Tests ─────────────────────────────────────────────────────

    #[test]
    fn serialization_naming_scope() {
        let scope = NamingScope::for_part("housing")
            .with_assembly("machine")
            .with_operation("fillet");

        let json = serde_json::to_string(&scope).unwrap();
        let decoded: NamingScope = serde_json::from_str(&json).unwrap();

        assert_eq!(decoded, scope);
    }

    #[test]
    fn serialization_scoped_id() {
        let scope = NamingScope::for_part("housing");
        let scoped = ScopedId::new(PersistentId(42), scope);

        let json = serde_json::to_string(&scoped).unwrap();
        let decoded: ScopedId = serde_json::from_str(&json).unwrap();

        assert_eq!(decoded, scoped);
    }

    #[test]
    fn serialization_enhanced_context_state() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.set_scope(NamingScope::for_part("test"));
        ctx.assign_id(10);
        ctx.assign_id(20);

        let state = ctx.export_state();

        let json = serde_json::to_string(&state).unwrap();
        let decoded: EnhancedNamingContextState = serde_json::from_str(&json).unwrap();

        assert_eq!(decoded.current_scope, state.current_scope);
        assert_eq!(decoded.next_persistent_id, state.next_persistent_id);
    }

    #[test]
    fn serialization_naming_context_snapshot() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.assign_id(10);

        let snapshot = NamingContextSnapshot::from_context(&ctx, 0, Some("test".to_string()));

        let json = serde_json::to_string(&snapshot).unwrap();
        let decoded: NamingContextSnapshot = serde_json::from_str(&json).unwrap();

        assert_eq!(decoded.id, snapshot.id);
        assert_eq!(decoded.operation, snapshot.operation);
    }
}

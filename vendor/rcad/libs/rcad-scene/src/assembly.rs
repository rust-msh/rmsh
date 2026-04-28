//! 装配体 / 实例化场景树
//!
//! # 设计原则
//!
//! - **`Arc<BRep>` 共享**：同一几何体可被多个 [`AssemblyNode`] 引用，不复制顶点数据。
//! - **`DAffine3` 变换叠加**：每个节点存储相对于父节点的仿射变换，查询时沿路径累积。
//! - **两种展开方式**：
//!   - [`Assembly::flatten`] — 返回 `(Arc<BRep>, 世界变换)` 列表，惰性实例化。
//!   - [`Assembly::to_brep`] — 合并为单一 BRep（调用 `BRep::transformed` + `append_brep`）。

use std::sync::Arc;
use std::collections::BTreeMap;

use glam::{DAffine3, DVec3};
use rcad_kernel::BRep;
use serde::{Deserialize, Serialize};

use crate::append_brep;

/// Semantic metadata carried by assembly/document nodes.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct AssemblyMetadata {
    /// Optional display label.
    pub display_name: Option<String>,
    /// Optional layer tag.
    pub layer: Option<String>,
    /// Optional material tag.
    pub material: Option<String>,
    /// Free-form key-value attributes.
    pub attributes: BTreeMap<String, String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// 核心类型
// ─────────────────────────────────────────────────────────────────────────────

/// 装配树中的单个节点。
///
/// 节点可以是叶节点（一个具体的 [`BRep`]）或子装配（包含若干子节点）。
/// `transform` 是相对于父节点的仿射变换，默认为恒等变换。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyNode {
    /// 节点唯一 ID（由所属 [`Assembly`] 分配）。
    pub id: u64,
    /// 人类可读名称。
    pub name: String,
    /// 相对于父节点的仿射变换（位置、旋转、缩放）。
    pub transform: DAffine3,
    /// 节点内容：叶节点几何或子装配列表。
    pub content: NodeContent,
    /// Semantic metadata for this node.
    #[serde(default)]
    pub metadata: AssemblyMetadata,
}

/// [`AssemblyNode`] 的内容类型。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum NodeContent {
    /// 叶节点：持有一个共享 BRep。
    Leaf(Arc<BRep>),
    /// 子装配：包含若干子节点（可以是叶节点或更深的子装配）。
    Assembly(Vec<AssemblyNode>),
}

/// 装配体根结构，持有零或多个顶层 [`AssemblyNode`]。
///
/// # 示例
///
/// ```
/// # use std::sync::Arc;
/// # use glam::{DAffine3, DVec3};
/// # use rcad_kernel::BRep;
/// use rcad_scene::assembly::Assembly;
///
/// let box_brep = Arc::new(BRep::new());
/// let mut asm = Assembly::new("my_assembly");
/// asm.add_part("part_a", Arc::clone(&box_brep));
/// asm.add_part_with_transform(
///     "part_b",
///     Arc::clone(&box_brep),
///     DAffine3::from_translation(DVec3::new(5.0, 0.0, 0.0)),
/// );
/// let flat = asm.flatten();
/// assert_eq!(flat.len(), 2);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Assembly {
    /// 装配体名称。
    pub name: String,
    /// 顶层节点列表（可以是叶节点或子装配）。
    pub roots: Vec<AssemblyNode>,
    /// Metadata for the whole assembly document.
    #[serde(default)]
    pub metadata: AssemblyMetadata,
    /// 下一个要分配的节点 ID。
    next_id: u64,
}

// ─────────────────────────────────────────────────────────────────────────────
// AssemblyNode 实现
// ─────────────────────────────────────────────────────────────────────────────

impl AssemblyNode {
    /// 创建叶节点（恒等变换）。
    pub fn new_leaf(id: u64, name: impl Into<String>, brep: Arc<BRep>) -> Self {
        Self {
            id,
            name: name.into(),
            transform: DAffine3::IDENTITY,
            content: NodeContent::Leaf(brep),
            metadata: AssemblyMetadata::default(),
        }
    }

    /// 创建带变换的叶节点。
    pub fn new_leaf_with_transform(
        id: u64,
        name: impl Into<String>,
        brep: Arc<BRep>,
        transform: DAffine3,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            transform,
            content: NodeContent::Leaf(brep),
            metadata: AssemblyMetadata::default(),
        }
    }

    /// 创建子装配节点（恒等变换，初始无子节点）。
    pub fn new_assembly(id: u64, name: impl Into<String>) -> Self {
        Self {
            id,
            name: name.into(),
            transform: DAffine3::IDENTITY,
            content: NodeContent::Assembly(Vec::new()),
            metadata: AssemblyMetadata::default(),
        }
    }

    /// 返回该节点下所有叶节点，累积 `parent_xform` 后写入 `out`。
    fn flatten_into(&self, parent_xform: DAffine3, out: &mut Vec<(Arc<BRep>, DAffine3)>) {
        let world = parent_xform * self.transform;
        match &self.content {
            NodeContent::Leaf(brep) => {
                out.push((Arc::clone(brep), world));
            }
            NodeContent::Assembly(children) => {
                for child in children {
                    child.flatten_into(world, out);
                }
            }
        }
    }

    /// 将该节点合并入目标 BRep（递归，累积父变换）。
    fn merge_into(&self, parent_xform: DAffine3, dst: &mut BRep) {
        let world = parent_xform * self.transform;
        match &self.content {
            NodeContent::Leaf(brep) => {
                let materialized = brep.transformed(world);
                append_brep(dst, materialized);
            }
            NodeContent::Assembly(children) => {
                for child in children {
                    child.merge_into(world, dst);
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Assembly 实现
// ─────────────────────────────────────────────────────────────────────────────

impl Assembly {
    /// 创建空装配体。
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            roots: Vec::new(),
            metadata: AssemblyMetadata::default(),
            next_id: 1,
        }
    }

    /// Attach or overwrite a document-level attribute.
    pub fn set_attribute(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.metadata.attributes.insert(key.into(), value.into());
    }

    /// Set assembly-level layer tag.
    pub fn set_layer(&mut self, layer: impl Into<String>) {
        self.metadata.layer = Some(layer.into());
    }

    /// Set assembly-level material tag.
    pub fn set_material(&mut self, material: impl Into<String>) {
        self.metadata.material = Some(material.into());
    }

    // ── 内部 ID 分配 ─────────────────────────────────────────────────────────

    fn alloc_id(&mut self) -> u64 {
        let id = self.next_id;
        self.next_id += 1;
        id
    }

    // ── 添加节点 ─────────────────────────────────────────────────────────────

    /// 添加顶层叶节点（恒等变换）。返回新节点的 ID。
    pub fn add_part(&mut self, name: impl Into<String>, brep: Arc<BRep>) -> u64 {
        let id = self.alloc_id();
        self.roots.push(AssemblyNode::new_leaf(id, name, brep));
        id
    }

    /// 添加顶层叶节点（带变换）。返回新节点的 ID。
    pub fn add_part_with_transform(
        &mut self,
        name: impl Into<String>,
        brep: Arc<BRep>,
        transform: DAffine3,
    ) -> u64 {
        let id = self.alloc_id();
        self.roots
            .push(AssemblyNode::new_leaf_with_transform(id, name, brep, transform));
        id
    }

    /// 添加一个已构建好的顶层节点（叶节点或子装配均可）。返回其 ID。
    pub fn add_node(&mut self, mut node: AssemblyNode) -> u64 {
        // 如果 id=0 表示未分配，则分配一个
        if node.id == 0 {
            node.id = self.alloc_id();
        } else {
            // 确保 next_id 不冲突
            if node.id >= self.next_id {
                self.next_id = node.id + 1;
            }
        }
        let id = node.id;
        self.roots.push(node);
        id
    }

    // ── 查询 ─────────────────────────────────────────────────────────────────

    /// 将整棵装配树展开为 `(共享 BRep, 世界变换)` 列表。
    ///
    /// 列表中每一项对应一个叶节点实例。变换已累积父链路，可直接用于渲染或坐标查询，
    /// 不复制 BRep 顶点数据。
    pub fn flatten(&self) -> Vec<(Arc<BRep>, DAffine3)> {
        let mut out = Vec::new();
        for root in &self.roots {
            root.flatten_into(DAffine3::IDENTITY, &mut out);
        }
        out
    }

    /// 将整棵装配树合并为单一 [`BRep`]。
    ///
    /// 对每个叶节点调用 [`BRep::transformed`]（只读副本），然后用 [`append_brep`] 拼接。
    /// 适用于需要单体 BRep 的场合（布尔运算、STEP 单体导出等）。
    pub fn to_brep(&self) -> BRep {
        let mut result = BRep::new();
        for root in &self.roots {
            root.merge_into(DAffine3::IDENTITY, &mut result);
        }
        result
    }

    /// 返回装配树中叶节点总数（即实例数）。
    pub fn instance_count(&self) -> usize {
        fn count(node: &AssemblyNode) -> usize {
            match &node.content {
                NodeContent::Leaf(_) => 1,
                NodeContent::Assembly(children) => children.iter().map(count).sum(),
            }
        }
        self.roots.iter().map(count).sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 从 AssemblyComponent (rcad-step) 构建 Assembly 的辅助
// ─────────────────────────────────────────────────────────────────────────────

/// 从一组扁平的 `(name, brep, transform)` 三元组快速构建 [`Assembly`]。
///
/// 这是将 `rcad_step::read_assembly` 返回值转换为场景树的推荐方式。
pub fn assembly_from_parts(
    name: impl Into<String>,
    parts: impl IntoIterator<Item = (String, BRep, DAffine3)>,
) -> Assembly {
    let mut asm = Assembly::new(name);
    for (part_name, brep, transform) in parts {
        asm.add_part_with_transform(part_name, Arc::new(brep), transform);
    }
    asm
}

// ─────────────────────────────────────────────────────────────────────────────
// DVec3 便利构造（减少调用侧 boilerplate）
// ─────────────────────────────────────────────────────────────────────────────

impl Assembly {
    /// 添加顶层叶节点，仅指定平移（旋转/缩放保持恒等）。
    pub fn add_part_at(
        &mut self,
        name: impl Into<String>,
        brep: Arc<BRep>,
        translation: DVec3,
    ) -> u64 {
        self.add_part_with_transform(name, brep, DAffine3::from_translation(translation))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assembly_metadata_setters_work() {
        let mut asm = Assembly::new("m");
        asm.set_layer("L1");
        asm.set_material("Aluminum");
        asm.set_attribute("owner", "team-a");

        assert_eq!(asm.metadata.layer.as_deref(), Some("L1"));
        assert_eq!(asm.metadata.material.as_deref(), Some("Aluminum"));
        assert_eq!(asm.metadata.attributes.get("owner").map(String::as_str), Some("team-a"));
    }

    #[test]
    fn node_metadata_defaults_empty() {
        let node = AssemblyNode::new_assembly(1, "sub");
        assert!(node.metadata.layer.is_none());
        assert!(node.metadata.attributes.is_empty());
    }
}

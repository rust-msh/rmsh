use glam::DAffine3;
use serde::{Deserialize, Serialize};

use crate::BRep;
use crate::appearance::Color;

/// 装配体组件引用：可以是一个独立 BRep 实体或嵌套的子装配体。
///
/// 类比 OCCT `XCAFDoc_ShapeTool` 的 shape reference。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ShapeRef {
    /// 叶节点：一个几何实体。
    Brep(BRep),
    /// 非叶节点：嵌套子装配体（box 避免递归类型无限大小）。
    Assembly(Box<Assembly>),
}

/// 装配体中的一个组件实例。
///
/// 同一个 `ShapeRef` 可被多个 `Component` 引用（实例化），
/// 每个实例有独立的变换和颜色覆盖。
///
/// 类比 OCCT `NEXT_ASSEMBLY_USAGE_OCCURENCE` + `XCAFDoc_Location`。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Component {
    /// 组件名称（用于选择、导出标注）。
    pub name: String,
    /// 几何内容（BRep 实体或子装配体）。
    pub shape: ShapeRef,
    /// 相对于父装配体坐标系的变换（位置 + 旋转 + 缩放）。
    pub transform: SerializableAffine3,
    /// 可选颜色覆盖（覆盖 shape 自身的颜色）。
    pub color: Option<Color>,
    /// 是否可见（用于渲染过滤）。
    pub visible: bool,
}

impl Component {
    /// 创建一个位于原点、无旋转的 BRep 组件。
    pub fn from_brep(name: impl Into<String>, brep: BRep) -> Self {
        Self {
            name: name.into(),
            shape: ShapeRef::Brep(brep),
            transform: SerializableAffine3::identity(),
            color: None,
            visible: true,
        }
    }

    /// 创建一个位于原点、无旋转的子装配体组件。
    pub fn from_assembly(name: impl Into<String>, asm: Assembly) -> Self {
        Self {
            name: name.into(),
            shape: ShapeRef::Assembly(Box::new(asm)),
            transform: SerializableAffine3::identity(),
            color: None,
            visible: true,
        }
    }

    /// 设置变换矩阵（Builder 风格）。
    pub fn with_transform(mut self, transform: DAffine3) -> Self {
        self.transform = SerializableAffine3(transform);
        self
    }

    /// 设置颜色覆盖（Builder 风格）。
    pub fn with_color(mut self, color: Color) -> Self {
        self.color = Some(color);
        self
    }

    /// 获取此组件的世界变换矩阵。
    pub fn affine(&self) -> DAffine3 {
        self.transform.0
    }
}

/// 层级装配体。
///
/// 可以包含多个 `Component`，每个组件有独立的变换。
/// 支持嵌套（装配体内嵌套子装配体）。
///
/// 类比 OCCT `XCAFDoc_ShapeTool` 管理的 shape 层级。
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Assembly {
    /// 装配体名称。
    pub name: String,
    /// 子组件列表（有序，用于渲染和遍历）。
    pub components: Vec<Component>,
}

impl Assembly {
    /// 创建一个空装配体。
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            components: Vec::new(),
        }
    }

    /// 添加一个组件。
    pub fn add(&mut self, component: Component) -> &mut Self {
        self.components.push(component);
        self
    }

    /// 返回此装配体中所有 BRep 叶节点（展平，带世界变换矩阵）。
    ///
    /// 用于渲染、导出、碰撞检测等需要遍历所有几何实体的场景。
    pub fn flatten(&self) -> Vec<FlatComponent> {
        let mut result = Vec::new();
        self.flatten_recursive(DAffine3::IDENTITY, &mut result);
        result
    }

    fn flatten_recursive(&self, parent_transform: DAffine3, out: &mut Vec<FlatComponent>) {
        for component in &self.components {
            if !component.visible {
                continue;
            }
            let world_transform = parent_transform * component.affine();
            match &component.shape {
                ShapeRef::Brep(brep) => {
                    out.push(FlatComponent {
                        name: component.name.clone(),
                        brep: brep.clone(),
                        world_transform,
                        color: component.color.clone(),
                    });
                }
                ShapeRef::Assembly(sub_asm) => {
                    sub_asm.flatten_recursive(world_transform, out);
                }
            }
        }
    }

    /// 递归查找指定名称的组件（深度优先）。
    pub fn find_component(&self, name: &str) -> Option<&Component> {
        for component in &self.components {
            if component.name == name {
                return Some(component);
            }
            if let ShapeRef::Assembly(sub) = &component.shape {
                if let Some(found) = sub.find_component(name) {
                    return Some(found);
                }
            }
        }
        None
    }

    /// 返回此装配体中直接子组件数量。
    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    /// 递归计算总叶节点（BRep）数量。
    pub fn leaf_count(&self) -> usize {
        self.components.iter().map(|c| match &c.shape {
            ShapeRef::Brep(_) => 1,
            ShapeRef::Assembly(sub) => sub.leaf_count(),
        }).sum()
    }
}

/// 展平后的叶节点（单个 BRep + 世界变换）。
#[derive(Debug, Clone)]
pub struct FlatComponent {
    pub name: String,
    pub brep: BRep,
    pub world_transform: DAffine3,
    pub color: Option<Color>,
}

/// `DAffine3` 的可序列化包装（glam 的 DAffine3 不直接支持 serde）。
#[derive(Debug, Clone, Copy)]
pub struct SerializableAffine3(pub DAffine3);

impl SerializableAffine3 {
    pub fn identity() -> Self {
        Self(DAffine3::IDENTITY)
    }
}

impl Serialize for SerializableAffine3 {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        let m = self.0.matrix3;
        let t = self.0.translation;
        // Serialize as 12 f64: [col0 x3, col1 x3, col2 x3, translation x3]
        let mut seq = serializer.serialize_seq(Some(12))?;
        for v in [m.x_axis.x, m.x_axis.y, m.x_axis.z,
                   m.y_axis.x, m.y_axis.y, m.y_axis.z,
                   m.z_axis.x, m.z_axis.y, m.z_axis.z,
                   t.x, t.y, t.z] {
            seq.serialize_element(&v)?;
        }
        seq.end()
    }
}

impl<'de> Deserialize<'de> for SerializableAffine3 {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let vals: Vec<f64> = Vec::deserialize(deserializer)?;
        if vals.len() != 12 {
            return Err(serde::de::Error::custom("expected 12 f64 for DAffine3"));
        }
        use glam::{DMat3, DVec3};
        let mat = DMat3::from_cols(
            DVec3::new(vals[0], vals[1], vals[2]),
            DVec3::new(vals[3], vals[4], vals[5]),
            DVec3::new(vals[6], vals[7], vals[8]),
        );
        let trans = DVec3::new(vals[9], vals[10], vals[11]);
        Ok(Self(DAffine3::from_mat3_translation(mat, trans)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BRep, PrimitiveSolid};

    #[test]
    fn assembly_flatten() {
        let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
        let box2 = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 0.5 });

        let mut asm = Assembly::new("root");
        asm.add(Component::from_brep("box", box1));
        asm.add(Component::from_brep("sphere", box2));

        let flat = asm.flatten();
        assert_eq!(flat.len(), 2);
        assert_eq!(flat[0].name, "box");
        assert_eq!(flat[1].name, "sphere");
    }

    #[test]
    fn assembly_nested() {
        let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });

        let mut sub = Assembly::new("sub");
        sub.add(Component::from_brep("part", box1));

        let mut root = Assembly::new("root");
        root.add(Component::from_assembly("sub_asm", sub));

        assert_eq!(root.leaf_count(), 1);
        let flat = root.flatten();
        assert_eq!(flat.len(), 1);
        assert_eq!(flat[0].name, "part");
    }

    #[test]
    fn assembly_find_component() {
        let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
        let mut asm = Assembly::new("root");
        asm.add(Component::from_brep("target", box1));

        assert!(asm.find_component("target").is_some());
        assert!(asm.find_component("missing").is_none());
    }
}

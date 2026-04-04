//! Geometry engine — bridges EMStudio's operation history with rcad's B-Rep kernel.
//!
//! The engine maintains an in-memory map of named objects to `BRep` instances.
//! It replays `GeometryOperation` steps by calling rcad APIs (primitives, booleans,
//! transforms, sweeps) and produces `GeoObject` snapshots for serialization.

use std::collections::HashMap;

use glam::{DAffine3, DMat3, DVec3};
use rcad_algorithms::{boolean_op, BooleanOpType};
use rcad_kernel::BRep;
use rcad_modeling::{box_brep, sphere_brep, cylinder_brep, cone_brep, torus_brep};
use thiserror::Error;

use crate::geometry::{BoundingBox, GeoObject, Geometry, GeometryOperation, ObjectAttributes, OperationCommand};

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

#[derive(Debug, Error)]
pub enum GeometryError {
    #[error("object not found: {0}")]
    ObjectNotFound(String),
    #[error("invalid parameters for {command}: {reason}")]
    InvalidParameters { command: String, reason: String },
    #[error("boolean operation failed: {0}")]
    BooleanFailed(String),
    #[error("build error: {0}")]
    BuildError(String),
    #[error("not implemented: {0}")]
    NotImplemented(String),
}

// ---------------------------------------------------------------------------
// Parameter extraction helpers
// ---------------------------------------------------------------------------

fn get_f64(params: &serde_json::Value, key: &str) -> Result<f64, GeometryError> {
    params
        .get(key)
        .and_then(|v| v.as_f64())
        .ok_or_else(|| GeometryError::InvalidParameters {
            command: String::new(),
            reason: format!("missing or invalid '{key}'"),
        })
}

fn get_vec3(params: &serde_json::Value, key: &str) -> Result<DVec3, GeometryError> {
    let arr = params
        .get(key)
        .and_then(|v| v.as_array())
        .ok_or_else(|| GeometryError::InvalidParameters {
            command: String::new(),
            reason: format!("missing or invalid '{key}'"),
        })?;
    if arr.len() < 3 {
        return Err(GeometryError::InvalidParameters {
            command: String::new(),
            reason: format!("'{key}' must be [x, y, z]"),
        });
    }
    Ok(DVec3::new(
        arr[0].as_f64().unwrap_or(0.0),
        arr[1].as_f64().unwrap_or(0.0),
        arr[2].as_f64().unwrap_or(0.0),
    ))
}

fn get_str<'a>(params: &'a serde_json::Value, key: &str) -> Result<&'a str, GeometryError> {
    params
        .get(key)
        .and_then(|v| v.as_str())
        .ok_or_else(|| GeometryError::InvalidParameters {
            command: String::new(),
            reason: format!("missing or invalid '{key}'"),
        })
}

fn get_vec3_or(params: &serde_json::Value, key: &str, default: DVec3) -> DVec3 {
    get_vec3(params, key).unwrap_or(default)
}

// ---------------------------------------------------------------------------
// GeometryEngine
// ---------------------------------------------------------------------------

/// Runtime geometry state — maintains named B-Rep objects derived from
/// the operation history.
pub struct GeometryEngine {
    /// Name → BRep mapping (runtime, not serialized).
    objects: HashMap<String, BRep>,
    /// Auto-incrementing object ID counter.
    next_id: u64,
}

impl GeometryEngine {
    pub fn new() -> Self {
        Self {
            objects: HashMap::new(),
            next_id: 1,
        }
    }

    /// Get a BRep by object name.
    pub fn get_brep(&self, name: &str) -> Option<&BRep> {
        self.objects.get(name)
    }

    /// Get all named BReps.
    pub fn all_breps(&self) -> &HashMap<String, BRep> {
        &self.objects
    }

    /// Rebuild all geometry by replaying an operation list from scratch.
    /// Returns the resulting `GeoObject` snapshot list.
    pub fn rebuild(
        &mut self,
        operations: &[GeometryOperation],
    ) -> Result<Vec<GeoObject>, GeometryError> {
        self.objects.clear();
        self.next_id = 1;

        let mut snapshots = Vec::new();

        for op in operations {
            let result = self.execute(op)?;
            if let Some(obj) = result {
                snapshots.push(obj);
            }
        }

        Ok(snapshots)
    }

    /// Execute a single geometry operation, returning a `GeoObject` snapshot
    /// if the operation creates or modifies a named object.
    pub fn execute(
        &mut self,
        op: &GeometryOperation,
    ) -> Result<Option<GeoObject>, GeometryError> {
        let params = &op.parameters;
        let attrs = op.attributes.as_ref().cloned().unwrap_or_default();

        match op.command {
            // -- Primitives --
            OperationCommand::CreateBox => {
                let origin = get_vec3_or(params, "origin", DVec3::ZERO);
                let size = get_vec3(params, "size").unwrap_or(DVec3::new(
                    get_f64(params, "width").unwrap_or(1.0),
                    get_f64(params, "height").unwrap_or(1.0),
                    get_f64(params, "depth").unwrap_or(1.0),
                ));
                let brep = box_brep(
                    origin,
                    DVec3::X,
                    DVec3::Y,
                    size.x,
                    size.y,
                    size.z,
                )
                .map_err(|e| GeometryError::BuildError(e.to_string()))?;

                let name = self.result_name(op, "Box");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::CreateCylinder => {
                let center = get_vec3_or(params, "center", DVec3::ZERO);
                let radius = get_f64(params, "radius")?;
                let height = get_f64(params, "height")?;
                let brep = cylinder_brep(
                    center,
                    DVec3::Z,
                    DVec3::X,
                    radius,
                    height,
                )
                .map_err(|e| GeometryError::BuildError(e.to_string()))?;

                let name = self.result_name(op, "Cylinder");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::CreateSphere => {
                let center = get_vec3_or(params, "center", DVec3::ZERO);
                let radius = get_f64(params, "radius")?;
                let brep = sphere_brep(center, radius)
                    .map_err(|e| GeometryError::BuildError(e.to_string()))?;

                let name = self.result_name(op, "Sphere");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::CreateCone => {
                let center = get_vec3_or(params, "center", DVec3::ZERO);
                let radius = get_f64(params, "radius").or_else(|_| get_f64(params, "base_radius"))?;
                let height = get_f64(params, "height")?;
                let brep = cone_brep(
                    center,
                    DVec3::Z,
                    DVec3::X,
                    radius,
                    height,
                )
                .map_err(|e| GeometryError::BuildError(e.to_string()))?;

                let name = self.result_name(op, "Cone");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::CreateTorus => {
                let center = get_vec3_or(params, "center", DVec3::ZERO);
                let major_radius = get_f64(params, "major_radius")?;
                let minor_radius = get_f64(params, "minor_radius")?;
                let brep = torus_brep(
                    center,
                    DVec3::Z,
                    DVec3::X,
                    major_radius,
                    minor_radius,
                )
                .map_err(|e| GeometryError::BuildError(e.to_string()))?;

                let name = self.result_name(op, "Torus");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            // -- Boolean Operations --
            OperationCommand::Unite => {
                self.execute_boolean(BooleanOpType::Union, params, op.step, &attrs)
            }
            OperationCommand::Subtract => {
                self.execute_boolean(BooleanOpType::Difference, params, op.step, &attrs)
            }
            OperationCommand::Intersect => {
                self.execute_boolean(BooleanOpType::Intersection, params, op.step, &attrs)
            }

            // -- Transforms --
            OperationCommand::Move => {
                let target = get_str(params, "target")?;
                let vector = get_vec3(params, "vector")?;
                let brep = self.take_brep(target)?;
                let mut brep = brep;
                brep.apply_transform(DAffine3::from_translation(vector));
                let name = target.to_string();
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::Rotate => {
                let target = get_str(params, "target")?;
                let axis = get_vec3_or(params, "axis", DVec3::Z);
                let angle_deg = get_f64(params, "angle_deg")?;
                let origin = get_vec3_or(params, "origin", DVec3::ZERO);
                let angle_rad = angle_deg.to_radians();
                let brep = self.take_brep(target)?;
                let mut brep = brep;
                // Translate to origin, rotate, translate back
                brep.apply_transform(DAffine3::from_translation(-origin));
                let rotation = DMat3::from_axis_angle(axis.normalize(), angle_rad);
                brep.apply_transform(DAffine3::from_mat3(rotation));
                brep.apply_transform(DAffine3::from_translation(origin));
                let name = target.to_string();
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::Scale => {
                let target = get_str(params, "target")?;
                let factor = get_f64(params, "factor").unwrap_or(1.0);
                let origin = get_vec3_or(params, "origin", DVec3::ZERO);
                let brep = self.take_brep(target)?;
                let mut brep = brep;
                brep.apply_transform(DAffine3::from_translation(-origin));
                brep.apply_transform(DAffine3::from_scale(DVec3::splat(factor)));
                brep.apply_transform(DAffine3::from_translation(origin));
                let name = target.to_string();
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::Mirror => {
                let target = get_str(params, "target")?;
                let normal = get_vec3_or(params, "normal", DVec3::X).normalize();
                let origin = get_vec3_or(params, "origin", DVec3::ZERO);
                let brep = self.take_brep(target)?;
                let mut brep = brep;
                // Reflection matrix: I - 2*n*nT
                let reflection = DMat3::from_cols(
                    DVec3::X - 2.0 * normal.x * normal,
                    DVec3::Y - 2.0 * normal.y * normal,
                    DVec3::Z - 2.0 * normal.z * normal,
                );
                brep.apply_transform(DAffine3::from_translation(-origin));
                brep.apply_transform(DAffine3::from_mat3(reflection));
                brep.apply_transform(DAffine3::from_translation(origin));
                let name = target.to_string();
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            // -- Property-only operations (no BRep change) --
            OperationCommand::SetMaterial
            | OperationCommand::SetColor
            | OperationCommand::SetSolveInside
            | OperationCommand::SetGroup => {
                // These only affect GeoObject metadata, not the BRep.
                // The attrs are stored in the operation and applied during snapshot.
                let target = get_str(params, "target")?;
                if !self.objects.contains_key(target) {
                    return Err(GeometryError::ObjectNotFound(target.to_string()));
                }
                // Re-snapshot with updated attrs
                let brep = self.objects.get(target).unwrap();
                let obj = self.make_snapshot(target, brep, op.step, &attrs);
                Ok(Some(obj))
            }

            OperationCommand::Rename => {
                let old_name = get_str(params, "old_name")
                    .or_else(|_| get_str(params, "target"))?;
                let new_name = get_str(params, "new_name")
                    .or_else(|_| get_str(params, "name"))?;
                let brep = self.take_brep(old_name)?;
                let name = new_name.to_string();
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            // -- Not yet implemented --
            OperationCommand::SweepAlongVector
            | OperationCommand::SweepAroundAxis
            | OperationCommand::SweepAlongPath
            | OperationCommand::CreatePolyline
            | OperationCommand::CreateRectangle
            | OperationCommand::CreateCircle
            | OperationCommand::DuplicateAlongLine
            | OperationCommand::DuplicateAroundAxis
            | OperationCommand::Fillet
            | OperationCommand::Chamfer
            | OperationCommand::Section
            | OperationCommand::Import => {
                Err(GeometryError::NotImplemented(format!("{:?}", op.command)))
            }
        }
    }

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------

    fn result_name(&self, op: &GeometryOperation, prefix: &str) -> String {
        op.result_object
            .clone()
            .unwrap_or_else(|| format!("{}{}", prefix, self.next_id))
    }

    fn insert(
        &mut self,
        name: String,
        brep: BRep,
        step: u32,
        attrs: &ObjectAttributes,
    ) -> GeoObject {
        let obj = self.make_snapshot(&name, &brep, step, attrs);
        self.objects.insert(name, brep);
        self.next_id += 1;
        obj
    }

    fn make_snapshot(
        &self,
        name: &str,
        brep: &BRep,
        step: u32,
        attrs: &ObjectAttributes,
    ) -> GeoObject {
        let bbox = brep.bounding_box().map(|[mn, mx]| BoundingBox {
            min: mn.to_array(),
            max: mx.to_array(),
        });
        GeoObject {
            id: self.next_id,
            name: name.to_string(),
            derived_from_step: step,
            material: attrs.material.clone().unwrap_or_else(|| "vacuum".to_string()),
            solve_inside: attrs.solve_inside.unwrap_or(false),
            color: attrs.color.unwrap_or([128, 128, 128]),
            transparency: attrs.transparency.unwrap_or(0.0),
            group: attrs.group.clone(),
            bounding_box: bbox,
        }
    }

    fn take_brep(&mut self, name: &str) -> Result<BRep, GeometryError> {
        self.objects
            .remove(name)
            .ok_or_else(|| GeometryError::ObjectNotFound(name.to_string()))
    }

    fn execute_boolean(
        &mut self,
        op_type: BooleanOpType,
        params: &serde_json::Value,
        step: u32,
        attrs: &ObjectAttributes,
    ) -> Result<Option<GeoObject>, GeometryError> {
        let target = get_str(params, "target")?;
        let tool = get_str(params, "tool")?;

        let brep_a = self.take_brep(target)?;
        let brep_b = self.take_brep(tool)?;

        let result = boolean_op(op_type, &brep_a, &brep_b)
            .map_err(|e| GeometryError::BooleanFailed(e.to_string()))?;

        let name = target.to_string();
        Ok(Some(self.insert(name, result, step, attrs)))
    }
}

impl Default for GeometryEngine {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Geometry integration
// ---------------------------------------------------------------------------

impl Geometry {
    /// Rebuild all geometry from the operation history using rcad.
    /// Returns a `GeometryEngine` that holds the live BRep objects.
    pub fn rebuild_with_engine(&mut self) -> Result<GeometryEngine, GeometryError> {
        let mut engine = GeometryEngine::new();
        self.objects = engine.rebuild(&self.operations)?;
        Ok(engine)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Geometry, GeometryOperation, OperationCommand};
    use serde_json::json;

    fn make_op(step: u32, cmd: OperationCommand, params: serde_json::Value) -> GeometryOperation {
        GeometryOperation {
            step,
            command: cmd,
            result_object: None,
            parameters: params,
            attributes: None,
        }
    }

    fn make_named_op(
        step: u32,
        cmd: OperationCommand,
        name: &str,
        params: serde_json::Value,
    ) -> GeometryOperation {
        GeometryOperation {
            step,
            command: cmd,
            result_object: Some(name.to_string()),
            parameters: params,
            attributes: None,
        }
    }

    #[test]
    fn create_box() {
        let mut engine = GeometryEngine::new();
        let op = make_named_op(
            1,
            OperationCommand::CreateBox,
            "Box1",
            json!({"origin": [0,0,0], "size": [10, 5, 3]}),
        );
        let result = engine.execute(&op).unwrap();
        assert!(result.is_some());
        let obj = result.unwrap();
        assert_eq!(obj.name, "Box1");
        assert!(obj.bounding_box.is_some());
        let bbox = obj.bounding_box.unwrap();
        // Box from origin [0,0,0] with size [10,5,3]
        assert!(bbox.max[0] > 0.0);
        assert!(engine.get_brep("Box1").is_some());
    }

    #[test]
    fn create_cylinder() {
        let mut engine = GeometryEngine::new();
        let op = make_named_op(
            1,
            OperationCommand::CreateCylinder,
            "Cyl1",
            json!({"center": [0,0,0], "radius": 5.0, "height": 10.0}),
        );
        let result = engine.execute(&op).unwrap().unwrap();
        assert_eq!(result.name, "Cyl1");
        assert!(result.bounding_box.is_some());
    }

    #[test]
    fn create_sphere() {
        let mut engine = GeometryEngine::new();
        let op = make_named_op(
            1,
            OperationCommand::CreateSphere,
            "Sph1",
            json!({"center": [0,0,0], "radius": 3.0}),
        );
        let result = engine.execute(&op).unwrap().unwrap();
        assert_eq!(result.name, "Sph1");
    }

    #[test]
    fn move_transform() {
        let mut engine = GeometryEngine::new();
        // Create a box at origin
        engine
            .execute(&make_named_op(
                1,
                OperationCommand::CreateBox,
                "Box1",
                json!({"origin": [0,0,0], "size": [2, 2, 2]}),
            ))
            .unwrap();

        // Move it by [10, 0, 0]
        let move_op = make_op(
            2,
            OperationCommand::Move,
            json!({"target": "Box1", "vector": [10, 0, 0]}),
        );
        let result = engine.execute(&move_op).unwrap().unwrap();
        let bbox = result.bounding_box.unwrap();
        // After moving by 10 in X, min x should be ~10
        assert!(bbox.min[0] >= 9.9, "min x should be ~10, got {}", bbox.min[0]);
    }

    #[test]
    fn rotate_transform() {
        let mut engine = GeometryEngine::new();
        engine
            .execute(&make_named_op(
                1,
                OperationCommand::CreateBox,
                "Box1",
                json!({"origin": [0,0,0], "size": [10, 2, 2]}),
            ))
            .unwrap();

        let rotate_op = make_op(
            2,
            OperationCommand::Rotate,
            json!({"target": "Box1", "axis": [0,0,1], "origin": [0,0,0], "angle_deg": 90}),
        );
        let result = engine.execute(&rotate_op).unwrap().unwrap();
        assert!(result.bounding_box.is_some());
    }

    #[test]
    fn boolean_unite() {
        let mut engine = GeometryEngine::new();
        engine
            .execute(&make_named_op(
                1,
                OperationCommand::CreateBox,
                "Box1",
                json!({"origin": [0,0,0], "size": [10, 10, 10]}),
            ))
            .unwrap();
        engine
            .execute(&make_named_op(
                2,
                OperationCommand::CreateBox,
                "Box2",
                json!({"origin": [5,0,0], "size": [10, 10, 10]}),
            ))
            .unwrap();

        let unite_op = make_op(
            3,
            OperationCommand::Unite,
            json!({"target": "Box1", "tool": "Box2"}),
        );
        let result = engine.execute(&unite_op).unwrap().unwrap();
        assert_eq!(result.name, "Box1");
        // Tool (Box2) should be consumed
        assert!(engine.get_brep("Box2").is_none());
        // Target (Box1) should still exist with the union result
        assert!(engine.get_brep("Box1").is_some());
    }

    #[test]
    fn object_not_found() {
        let mut engine = GeometryEngine::new();
        let op = make_op(
            1,
            OperationCommand::Move,
            json!({"target": "NonExistent", "vector": [1,0,0]}),
        );
        let result = engine.execute(&op);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), GeometryError::ObjectNotFound(_)));
    }

    #[test]
    fn rebuild_from_operations() {
        let ops = vec![
            make_named_op(
                1,
                OperationCommand::CreateBox,
                "Box1",
                json!({"origin": [0,0,0], "size": [10, 10, 10]}),
            ),
            make_named_op(
                2,
                OperationCommand::CreateCylinder,
                "Cyl1",
                json!({"center": [5,5,0], "radius": 3.0, "height": 10.0}),
            ),
            make_op(
                3,
                OperationCommand::Move,
                json!({"target": "Box1", "vector": [20, 0, 0]}),
            ),
        ];

        let mut geom = Geometry {
            operations: ops,
            objects: Vec::new(),
        };
        let engine = geom.rebuild_with_engine().unwrap();
        assert_eq!(geom.objects.len(), 3); // Box1 created, Cyl1 created, Box1 moved (updated)
        assert!(engine.get_brep("Box1").is_some());
        assert!(engine.get_brep("Cyl1").is_some());
    }

    #[test]
    fn rename_object() {
        let mut engine = GeometryEngine::new();
        engine
            .execute(&make_named_op(
                1,
                OperationCommand::CreateBox,
                "Box1",
                json!({"origin": [0,0,0], "size": [5, 5, 5]}),
            ))
            .unwrap();

        let rename_op = make_op(
            2,
            OperationCommand::Rename,
            json!({"target": "Box1", "name": "MyBox"}),
        );
        engine.execute(&rename_op).unwrap();
        assert!(engine.get_brep("Box1").is_none());
        assert!(engine.get_brep("MyBox").is_some());
    }
}

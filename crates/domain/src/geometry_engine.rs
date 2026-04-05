//! Geometry engine — bridges EMStudio's operation history with rcad's B-Rep kernel.
//!
//! The engine maintains an in-memory map of named objects to `BRep` instances.
//! It replays `GeometryOperation` steps by calling rcad APIs (primitives, booleans,
//! transforms, sweeps) and produces `GeoObject` snapshots for serialization.

use std::collections::HashMap;

use glam::{DAffine3, DMat3, DVec2, DVec3};
use rcad_algorithms::{boolean_op, BooleanOpType};
use rcad_kernel::BRep;
use rcad_modeling::{box_brep, sphere_brep, cylinder_brep, cone_brep, torus_brep, extrude, revolve, sweep_pipe};
use thiserror::Error;

use crate::expression::{evaluate, parse_value_with_unit};
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
// Parameter extraction helpers (parametric — supports expressions)
// ---------------------------------------------------------------------------

/// Evaluate a parameter as f64. Supports:
/// - JSON number: `10.0`
/// - JSON string with unit: `"10mm"`, `"2.4GHz"`
/// - JSON string expression with variables: `"antenna_width * 2"`
fn eval_f64(
    params: &serde_json::Value,
    key: &str,
    vars: &HashMap<String, f64>,
) -> Result<f64, GeometryError> {
    let val = params.get(key).ok_or_else(|| GeometryError::InvalidParameters {
        command: String::new(),
        reason: format!("missing '{key}'"),
    })?;

    // Direct number
    if let Some(n) = val.as_f64() {
        return Ok(n);
    }
    // Integer stored as i64
    if let Some(n) = val.as_i64() {
        return Ok(n as f64);
    }
    // String: try value-with-unit, then expression
    if let Some(s) = val.as_str() {
        return parse_value_with_unit(s)
            .or_else(|_| evaluate(s, vars))
            .map_err(|e| GeometryError::InvalidParameters {
                command: String::new(),
                reason: format!("cannot evaluate '{key}': {e}"),
            });
    }

    Err(GeometryError::InvalidParameters {
        command: String::new(),
        reason: format!("'{key}' is not a number or evaluable string"),
    })
}

/// Evaluate a vec3 parameter. Supports:
/// - JSON array of numbers: `[1.0, 2.0, 3.0]`
/// - JSON array of strings: `["10mm", "antenna_width", "0"]`
fn eval_vec3(
    params: &serde_json::Value,
    key: &str,
    vars: &HashMap<String, f64>,
) -> Result<DVec3, GeometryError> {
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

    let resolve_component = |v: &serde_json::Value, idx: usize| -> Result<f64, GeometryError> {
        if let Some(n) = v.as_f64() {
            return Ok(n);
        }
        if let Some(n) = v.as_i64() {
            return Ok(n as f64);
        }
        if let Some(s) = v.as_str() {
            return parse_value_with_unit(s)
                .or_else(|_| evaluate(s, vars))
                .map_err(|e| GeometryError::InvalidParameters {
                    command: String::new(),
                    reason: format!("cannot evaluate '{key}[{idx}]': {e}"),
                });
        }
        Ok(0.0)
    };

    Ok(DVec3::new(
        resolve_component(&arr[0], 0)?,
        resolve_component(&arr[1], 1)?,
        resolve_component(&arr[2], 2)?,
    ))
}

fn eval_vec3_or(
    params: &serde_json::Value,
    key: &str,
    vars: &HashMap<String, f64>,
    default: DVec3,
) -> DVec3 {
    eval_vec3(params, key, vars).unwrap_or(default)
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
    /// `vars` contains resolved variable values for parametric evaluation.
    /// Returns the resulting `GeoObject` snapshot list.
    pub fn rebuild(
        &mut self,
        operations: &[GeometryOperation],
        vars: &HashMap<String, f64>,
    ) -> Result<Vec<GeoObject>, GeometryError> {
        self.objects.clear();
        self.next_id = 1;

        let mut snapshots = Vec::new();

        for op in operations {
            let result = self.execute(op, vars)?;
            if let Some(obj) = result {
                snapshots.push(obj);
            }
        }

        Ok(snapshots)
    }

    /// Execute a single geometry operation, returning a `GeoObject` snapshot
    /// if the operation creates or modifies a named object.
    /// `vars` contains resolved variable values for parametric evaluation.
    pub fn execute(
        &mut self,
        op: &GeometryOperation,
        vars: &HashMap<String, f64>,
    ) -> Result<Option<GeoObject>, GeometryError> {
        let params = &op.parameters;
        let attrs = op.attributes.as_ref().cloned().unwrap_or_default();

        match op.command {
            // -- Primitives --
            OperationCommand::CreateBox => {
                let origin = eval_vec3_or(params, "origin", vars, DVec3::ZERO);
                let size = eval_vec3(params, "size", vars).unwrap_or(DVec3::new(
                    eval_f64(params, "width", vars).unwrap_or(1.0),
                    eval_f64(params, "height", vars).unwrap_or(1.0),
                    eval_f64(params, "depth", vars).unwrap_or(1.0),
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
                let center = eval_vec3_or(params, "center", vars, DVec3::ZERO);
                let radius = eval_f64(params, "radius", vars)?;
                let height = eval_f64(params, "height", vars)?;
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
                let center = eval_vec3_or(params, "center", vars, DVec3::ZERO);
                let radius = eval_f64(params, "radius", vars)?;
                let brep = sphere_brep(center, radius)
                    .map_err(|e| GeometryError::BuildError(e.to_string()))?;

                let name = self.result_name(op, "Sphere");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::CreateCone => {
                let center = eval_vec3_or(params, "center", vars, DVec3::ZERO);
                let radius = eval_f64(params, "radius", vars).or_else(|_| eval_f64(params, "base_radius", vars))?;
                let height = eval_f64(params, "height", vars)?;
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
                let center = eval_vec3_or(params, "center", vars, DVec3::ZERO);
                let major_radius = eval_f64(params, "major_radius", vars)?;
                let minor_radius = eval_f64(params, "minor_radius", vars)?;
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
                let vector = eval_vec3(params, "vector", vars)?;
                let brep = self.take_brep(target)?;
                let mut brep = brep;
                brep.apply_transform(DAffine3::from_translation(vector));
                let name = target.to_string();
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::Rotate => {
                let target = get_str(params, "target")?;
                let axis = eval_vec3_or(params, "axis", vars, DVec3::Z);
                let angle_deg = eval_f64(params, "angle_deg", vars)?;
                let origin = eval_vec3_or(params, "origin", vars, DVec3::ZERO);
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
                let factor = eval_f64(params, "factor", vars).unwrap_or(1.0);
                let origin = eval_vec3_or(params, "origin", vars, DVec3::ZERO);
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
                let normal = eval_vec3_or(params, "normal", vars, DVec3::X).normalize();
                let origin = eval_vec3_or(params, "origin", vars, DVec3::ZERO);
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

            // -- Sweeps --
            OperationCommand::SweepAlongVector => {
                let target = get_str(params, "target")?;
                let direction = eval_vec3(params, "direction", vars)?;
                let distance = eval_f64(params, "distance", vars)?;
                let face_idx = params
                    .get("face_index")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(0) as usize;
                let profile = self.take_brep(target)?;
                let brep = extrude(&profile, face_idx, direction, distance)
                    .map_err(|e| GeometryError::BuildError(e.to_string()))?;
                let name = self.result_name(op, "Extrude");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::SweepAroundAxis => {
                let target = get_str(params, "target")?;
                let axis_origin = eval_vec3_or(params, "axis_origin", vars, DVec3::ZERO);
                let axis_direction = eval_vec3(params, "axis_direction", vars)?;
                let angle_deg = eval_f64(params, "angle", vars)?;
                let angle_rad = angle_deg.to_radians();
                let face_idx = params
                    .get("face_index")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(0) as usize;
                let profile = self.take_brep(target)?;
                let brep = revolve(&profile, face_idx, axis_origin, axis_direction, angle_rad)
                    .map_err(|e| GeometryError::BuildError(e.to_string()))?;
                let name = self.result_name(op, "Revolve");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::SweepAlongPath => {
                let profile_arr = params
                    .get("profile")
                    .and_then(|v| v.as_array())
                    .ok_or_else(|| GeometryError::InvalidParameters {
                        command: "SweepAlongPath".to_string(),
                        reason: "missing 'profile' array of [x, y] points".to_string(),
                    })?;
                let profile_2d: Vec<DVec2> = profile_arr
                    .iter()
                    .filter_map(|p| {
                        let a = p.as_array()?;
                        Some(DVec2::new(a.first()?.as_f64()?, a.get(1)?.as_f64()?))
                    })
                    .collect();

                let spine_arr = params
                    .get("spine")
                    .and_then(|v| v.as_array())
                    .ok_or_else(|| GeometryError::InvalidParameters {
                        command: "SweepAlongPath".to_string(),
                        reason: "missing 'spine' array of [x, y, z] points".to_string(),
                    })?;
                let spine: Vec<DVec3> = spine_arr
                    .iter()
                    .filter_map(|p| {
                        let a = p.as_array()?;
                        Some(DVec3::new(
                            a.first()?.as_f64()?,
                            a.get(1)?.as_f64()?,
                            a.get(2)?.as_f64()?,
                        ))
                    })
                    .collect();

                let brep = sweep_pipe(&profile_2d, &spine)
                    .map_err(|e| GeometryError::BuildError(e.to_string()))?;
                let name = self.result_name(op, "Sweep");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            OperationCommand::Import => {
                let file_path = get_str(params, "file_path")?;
                let brep = rcad_step::StepReader::read_file(file_path)
                    .map_err(|e| GeometryError::BuildError(e))?;
                let name = self.result_name(op, "Import");
                Ok(Some(self.insert(name, brep, op.step, &attrs)))
            }

            // -- Not yet implemented --
            OperationCommand::CreatePolyline
            | OperationCommand::CreateRectangle
            | OperationCommand::CreateCircle
            | OperationCommand::DuplicateAlongLine
            | OperationCommand::DuplicateAroundAxis
            | OperationCommand::Fillet
            | OperationCommand::Chamfer
            | OperationCommand::Section => {
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
    /// `vars` contains resolved variable values for parametric evaluation.
    /// Returns a `GeometryEngine` that holds the live BRep objects.
    pub fn rebuild_with_engine(
        &mut self,
        vars: &HashMap<String, f64>,
    ) -> Result<GeometryEngine, GeometryError> {
        let mut engine = GeometryEngine::new();
        self.objects = engine.rebuild(&self.operations, vars)?;
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
        let vars = HashMap::new();
        let op = make_named_op(
            1,
            OperationCommand::CreateBox,
            "Box1",
            json!({"origin": [0,0,0], "size": [10, 5, 3]}),
        );
        let result = engine.execute(&op, &vars).unwrap();
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
        let vars = HashMap::new();
        let op = make_named_op(
            1,
            OperationCommand::CreateCylinder,
            "Cyl1",
            json!({"center": [0,0,0], "radius": 5.0, "height": 10.0}),
        );
        let result = engine.execute(&op, &vars).unwrap().unwrap();
        assert_eq!(result.name, "Cyl1");
        assert!(result.bounding_box.is_some());
    }

    #[test]
    fn create_sphere() {
        let mut engine = GeometryEngine::new();
        let vars = HashMap::new();
        let op = make_named_op(
            1,
            OperationCommand::CreateSphere,
            "Sph1",
            json!({"center": [0,0,0], "radius": 3.0}),
        );
        let result = engine.execute(&op, &vars).unwrap().unwrap();
        assert_eq!(result.name, "Sph1");
    }

    #[test]
    fn move_transform() {
        let mut engine = GeometryEngine::new();
        let vars = HashMap::new();
        engine
            .execute(&make_named_op(
                1,
                OperationCommand::CreateBox,
                "Box1",
                json!({"origin": [0,0,0], "size": [2, 2, 2]}),
            ), &vars)
            .unwrap();

        let move_op = make_op(
            2,
            OperationCommand::Move,
            json!({"target": "Box1", "vector": [10, 0, 0]}),
        );
        let result = engine.execute(&move_op, &vars).unwrap().unwrap();
        let bbox = result.bounding_box.unwrap();
        assert!(bbox.min[0] >= 9.9, "min x should be ~10, got {}", bbox.min[0]);
    }

    #[test]
    fn rotate_transform() {
        let mut engine = GeometryEngine::new();
        let vars = HashMap::new();
        engine
            .execute(&make_named_op(
                1,
                OperationCommand::CreateBox,
                "Box1",
                json!({"origin": [0,0,0], "size": [10, 2, 2]}),
            ), &vars)
            .unwrap();

        let rotate_op = make_op(
            2,
            OperationCommand::Rotate,
            json!({"target": "Box1", "axis": [0,0,1], "origin": [0,0,0], "angle_deg": 90}),
        );
        let result = engine.execute(&rotate_op, &vars).unwrap().unwrap();
        assert!(result.bounding_box.is_some());
    }

    #[test]
    fn boolean_unite() {
        let mut engine = GeometryEngine::new();
        let vars = HashMap::new();
        engine
            .execute(&make_named_op(
                1,
                OperationCommand::CreateBox,
                "Box1",
                json!({"origin": [0,0,0], "size": [10, 10, 10]}),
            ), &vars)
            .unwrap();
        engine
            .execute(&make_named_op(
                2,
                OperationCommand::CreateBox,
                "Box2",
                json!({"origin": [5,0,0], "size": [10, 10, 10]}),
            ), &vars)
            .unwrap();

        let unite_op = make_op(
            3,
            OperationCommand::Unite,
            json!({"target": "Box1", "tool": "Box2"}),
        );
        let result = engine.execute(&unite_op, &vars).unwrap().unwrap();
        assert_eq!(result.name, "Box1");
        assert!(engine.get_brep("Box2").is_none());
        assert!(engine.get_brep("Box1").is_some());
    }

    #[test]
    fn object_not_found() {
        let mut engine = GeometryEngine::new();
        let vars = HashMap::new();
        let op = make_op(
            1,
            OperationCommand::Move,
            json!({"target": "NonExistent", "vector": [1,0,0]}),
        );
        let result = engine.execute(&op, &vars);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), GeometryError::ObjectNotFound(_)));
    }

    #[test]
    fn rebuild_from_operations() {
        let vars = HashMap::new();
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
        let engine = geom.rebuild_with_engine(&vars).unwrap();
        assert_eq!(geom.objects.len(), 3);
        assert!(engine.get_brep("Box1").is_some());
        assert!(engine.get_brep("Cyl1").is_some());
    }

    #[test]
    fn rename_object() {
        let mut engine = GeometryEngine::new();
        let vars = HashMap::new();
        engine
            .execute(&make_named_op(
                1,
                OperationCommand::CreateBox,
                "Box1",
                json!({"origin": [0,0,0], "size": [5, 5, 5]}),
            ), &vars)
            .unwrap();

        let rename_op = make_op(
            2,
            OperationCommand::Rename,
            json!({"target": "Box1", "name": "MyBox"}),
        );
        engine.execute(&rename_op, &vars).unwrap();
        assert!(engine.get_brep("Box1").is_none());
        assert!(engine.get_brep("MyBox").is_some());
    }

    #[test]
    fn parametric_box_with_string_dimensions() {
        let mut engine = GeometryEngine::new();
        let vars = HashMap::new();
        let op = make_named_op(
            1,
            OperationCommand::CreateBox,
            "Box1",
            json!({"origin": [0,0,0], "width": "10mm", "height": "5mm", "depth": "3mm"}),
        );
        let result = engine.execute(&op, &vars).unwrap().unwrap();
        assert_eq!(result.name, "Box1");
        assert!(result.bounding_box.is_some());
    }

    #[test]
    fn parametric_box_with_variable_reference() {
        let mut engine = GeometryEngine::new();
        let mut vars = HashMap::new();
        vars.insert("w".to_string(), 10.0);
        vars.insert("h".to_string(), 5.0);
        vars.insert("d".to_string(), 3.0);
        let op = make_named_op(
            1,
            OperationCommand::CreateBox,
            "Box1",
            json!({"origin": [0,0,0], "width": "w", "height": "h", "depth": "d"}),
        );
        let result = engine.execute(&op, &vars).unwrap().unwrap();
        assert_eq!(result.name, "Box1");
        assert!(result.bounding_box.is_some());
    }

    #[test]
    fn parametric_box_with_expression() {
        let mut engine = GeometryEngine::new();
        let mut vars = HashMap::new();
        vars.insert("base".to_string(), 5.0);
        let op = make_named_op(
            1,
            OperationCommand::CreateBox,
            "Box1",
            json!({"origin": [0,0,0], "width": "base * 2", "height": "base", "depth": "base / 2"}),
        );
        let result = engine.execute(&op, &vars).unwrap().unwrap();
        assert_eq!(result.name, "Box1");
        let bbox = result.bounding_box.unwrap();
        // width=10, height=5, depth=2.5
        assert!((bbox.max[0] - 10.0).abs() < 0.1);
        assert!((bbox.max[1] - 5.0).abs() < 0.1);
        assert!((bbox.max[2] - 2.5).abs() < 0.1);
    }
}

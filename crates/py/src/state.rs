//! Global runtime state for the gmsh-style Python API.
//!
//! The gmsh API is procedural with process-global state (`initialize`,
//! `model.occ.addBox`, ...), so the bindings keep one global registry of
//! named models behind a mutex. A [`Model`] carries the CAD entities
//! created through `model.occ.*` (real `rcad_kernel::BRep` solids built
//! by rcad-modeling primitives and the OCCT-aligned rcad-algo booleans),
//! the synchronized [`rmsh_model::Mesh`], physical groups and the option
//! overrides.

use std::collections::BTreeMap;
use std::sync::{LazyLock, Mutex};

use glam::DVec3;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::PyResult;

/// One CAD entity of a model.
pub enum OccShape {
    /// A real B-Rep solid or face (rcad kernel).
    BRep { dim: i32, brep: rcad_kernel::BRep },
    /// A point entity (construction geometry; meshing ignores it).
    Point(DVec3),
    /// A line segment entity (construction geometry).
    Line(DVec3, DVec3),
    /// A circular arc through start-center-end (construction geometry).
    Arc(DVec3, DVec3, DVec3),
    /// A spline through the given points (construction geometry).
    Spline(Vec<DVec3>),
    /// Curve / surface loop bookkeeping (lists of entity tags).
    Loop { dim: i32, members: Vec<i32> },
}

impl OccShape {
    pub fn dimension(&self) -> i32 {
        match self {
            OccShape::BRep { dim, .. } => *dim,
            OccShape::Point(_) => 0,
            OccShape::Line(..) | OccShape::Arc(..) | OccShape::Spline(_) => 1,
            OccShape::Loop { dim, .. } => *dim,
        }
    }
}

pub struct PhysicalGroup {
    pub dim: i32,
    pub tag: i32,
    pub entities: Vec<i32>,
    #[allow(dead_code)] // consumed through mesh.physical_names at write time
    pub name: Option<String>,
}

pub struct Model {
    #[allow(dead_code)] // the registry keys by name; the field stays for introspection
    pub name: String,
    pub shapes: BTreeMap<i32, OccShape>,
    pub next_tag: i32,
    /// The synchronized + generated mesh (Triangle3 surface, tets after
    /// `mesh.generate(3)`).
    pub mesh: rmsh_model::Mesh,
    /// Maps element id -> the CAD entity tag it was synchronized from, so
    /// physical groups can be applied at write time.
    pub element_entity: BTreeMap<u64, i32>,
    pub physicals: Vec<PhysicalGroup>,
    pub next_physical_tag: i32,
    /// Per-entity mesh size from `mesh.setSize` (v1: the minimum over the
    /// entries applies globally at generate time).
    pub mesh_size: Option<f64>,
}

impl Model {
    pub fn new(name: String) -> Self {
        Model {
            name,
            shapes: BTreeMap::new(),
            next_tag: 1,
            mesh: rmsh_model::Mesh::new(),
            element_entity: BTreeMap::new(),
            physicals: Vec::new(),
            next_physical_tag: 1,
            mesh_size: None,
        }
    }

    pub fn alloc_tag(&mut self) -> i32 {
        let tag = self.next_tag;
        self.next_tag += 1;
        tag
    }

    pub fn bbox(&self) -> ([f64; 3], [f64; 3]) {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for shape in self.shapes.values() {
            match shape {
                OccShape::BRep { brep, .. } => {
                    if let Some([bmin, bmax]) = brep.bounding_box() {
                        for d in 0..3 {
                            min[d] = min[d].min(bmin[d]);
                            max[d] = max[d].max(bmax[d]);
                        }
                    }
                }
                OccShape::Point(p) => {
                    for d in 0..3 {
                        min[d] = min[d].min(p[d]);
                        max[d] = max[d].max(p[d]);
                    }
                }
                OccShape::Line(a, b) | OccShape::Arc(a, _, b) => {
                    for p in [a, b] {
                        for d in 0..3 {
                            min[d] = min[d].min(p[d]);
                            max[d] = max[d].max(p[d]);
                        }
                    }
                }
                OccShape::Spline(pts) => {
                    for p in pts {
                        for d in 0..3 {
                            min[d] = min[d].min(p[d]);
                            max[d] = max[d].max(p[d]);
                        }
                    }
                }
                OccShape::Loop { .. } => {}
            }
        }
        (min, max)
    }
}

#[derive(Clone, Debug)]
pub enum OptionValue {
    Number(f64),
    String(String),
    Color((i32, i32, i32, i32)),
}

pub struct RuntimeState {
    pub initialized: bool,
    pub models: BTreeMap<String, Model>,
    pub current: Option<String>,
    pub logs: Vec<String>,
    pub logger_on: bool,
    pub options: BTreeMap<String, OptionValue>,
}

impl Default for RuntimeState {
    fn default() -> Self {
        RuntimeState {
            initialized: false,
            models: BTreeMap::new(),
            current: None,
            logs: Vec::new(),
            logger_on: false,
            options: BTreeMap::new(),
        }
    }
}

static STATE: LazyLock<Mutex<RuntimeState>> = LazyLock::new(|| Mutex::new(RuntimeState::default()));

pub fn lock() -> std::sync::MutexGuard<'static, RuntimeState> {
    STATE.lock().unwrap_or_else(|poisoned| poisoned.into_inner())
}

pub fn require_initialized(state: &RuntimeState) -> PyResult<()> {
    if state.initialized {
        Ok(())
    } else {
        Err(PyRuntimeError::new_err(
            "rmsh is not initialized; call initialize() first",
        ))
    }
}

/// Run `f` with the current model mutable.
pub fn with_current_model<T>(
    f: impl FnOnce(&mut Model) -> PyResult<T>,
) -> PyResult<T> {
    let mut state = lock();
    require_initialized(&state)?;
    let name = state
        .current
        .clone()
        .ok_or_else(|| PyRuntimeError::new_err("no current model; call model.add() first"))?;
    let model = state
        .models
        .get_mut(&name)
        .ok_or_else(|| PyRuntimeError::new_err(format!("current model {name:?} vanished")))?;
    f(model)
}

/// Run `f` with the current model shared.
pub fn with_current_model_ref<T>(
    f: impl FnOnce(&Model) -> PyResult<T>,
) -> PyResult<T> {
    let state = lock();
    require_initialized(&state)?;
    let name = state
        .current
        .clone()
        .ok_or_else(|| PyRuntimeError::new_err("no current model; call model.add() first"))?;
    let model = state
        .models
        .get(&name)
        .ok_or_else(|| PyRuntimeError::new_err(format!("current model {name:?} vanished")))?;
    f(model)
}

/// Run `f` with the option table shared and the current model mutable —
/// for generate paths that read options while meshing.
pub fn with_current_model_options<T>(
    f: impl FnOnce(&BTreeMap<String, OptionValue>, &mut Model) -> PyResult<T>,
) -> PyResult<T> {
    let mut state = lock();
    require_initialized(&state)?;
    let name = state
        .current
        .clone()
        .ok_or_else(|| PyRuntimeError::new_err("no current model; call model.add() first"))?;
    let options = state.options.clone();
    let model = state
        .models
        .get_mut(&name)
        .ok_or_else(|| PyRuntimeError::new_err(format!("current model {name:?} vanished")))?;
    f(&options, model)
}

pub fn option_number(options: &BTreeMap<String, OptionValue>, name: &str) -> Option<f64> {
    match options.get(name) {
        Some(OptionValue::Number(v)) => Some(*v),
        Some(_) => None,
        None => match default_option(name) {
            Some(OptionValue::Number(v)) => Some(v),
            _ => None,
        },
    }
}

/// The option names the bindings honor. `gmsh` rejects unknown names and so
/// do we — silently accepting them would fake compatibility.
pub const KNOWN_OPTIONS: &[(&str, OptionValue)] = &[
    ("Mesh.MeshSizeMax", OptionValue::Number(1e22)),
    ("Mesh.MeshSizeMin", OptionValue::Number(0.0)),
    ("Mesh.CharacteristicLengthMax", OptionValue::Number(1e22)),
    ("Mesh.CharacteristicLengthMin", OptionValue::Number(0.0)),
    ("Mesh.MshFileVersion", OptionValue::Number(2.2)),
    ("Mesh.Binary", OptionValue::Number(0.0)),
    ("Mesh.Algorithm", OptionValue::Number(6.0)),
    ("Mesh.Algorithm3D", OptionValue::Number(1.0)),
    ("General.Verbosity", OptionValue::Number(5.0)),
];

pub fn default_option(name: &str) -> Option<OptionValue> {
    KNOWN_OPTIONS
        .iter()
        .find(|(n, _)| *n == name)
        .map(|(_, v)| v.clone())
}

pub fn check_option_name(name: &str) -> PyResult<()> {
    if default_option(name).is_some() {
        Ok(())
    } else {
        Err(PyValueError::new_err(format!(
            "unknown option {name:?}; supported: {}",
            KNOWN_OPTIONS
                .iter()
                .map(|(n, _)| *n)
                .collect::<Vec<_>>()
                .join(", ")
        )))
    }
}

#[allow(dead_code)] // the logger surface lands with the quality-stat binding
pub fn log(state: &mut RuntimeState, message: impl Into<String>) {
    state.logs.push(message.into());
    if state.logs.len() > 4096 {
        state.logs.remove(0);
    }
}

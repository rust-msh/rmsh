use serde::{Deserialize, Serialize};

/// PCB / package layer stack definition.
///
/// Layers are ordered from bottom to top. The first layer is the bottom-most
/// (typically a ground plane or air half-space), and the last layer is the
/// top-most (typically air half-space).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerStack {
    pub name: String,
    /// Stack unit (default: "mm").
    #[serde(default = "default_stack_unit")]
    pub unit: String,
    /// Layers from bottom to top.
    pub layers: Vec<Layer>,
}

fn default_stack_unit() -> String {
    "mm".to_string()
}

/// A single layer in the PCB stackup.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Layer {
    pub name: String,
    pub layer_type: LayerType,
    /// Layer thickness in stack units.
    pub thickness: f64,
    /// Z offset from the global origin (bottom of this layer).
    #[serde(default)]
    pub z_offset: f64,
    /// Material name reference (must exist in Design.definitions.materials).
    #[serde(default)]
    pub material: Option<String>,
    /// Relative permittivity override (if no material assigned).
    #[serde(default)]
    pub eps_r: Option<f64>,
    /// Loss tangent override.
    #[serde(default)]
    pub loss_tangent: Option<f64>,
    /// Conductivity override [S/m] (for conductor layers).
    #[serde(default)]
    pub conductivity: Option<f64>,
}

/// Classification of a layer in the stackup.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum LayerType {
    /// Signal trace layer (copper).
    Signal,
    /// Ground / power plane layer.
    Ground,
    /// Dielectric / substrate layer (FR4, Rogers, etc.).
    Dielectric,
    /// Via / through-hole layer.
    Via,
    /// Solder mask layer.
    SolderMask,
    /// Silkscreen layer.
    Silkscreen,
    /// Background (air half-space) — typically the outermost layers.
    Background,
}

impl LayerStack {
    pub fn new(name: impl Into<String>, unit: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            unit: unit.into(),
            layers: Vec::new(),
        }
    }

    /// Add a layer and auto-compute its z_offset from the previous layer.
    pub fn push_layer(&mut self, mut layer: Layer) {
        let z = self.total_thickness();
        layer.z_offset = z;
        self.layers.push(layer);
    }

    /// Total thickness of all layers (bottom to top).
    pub fn total_thickness(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness).sum()
    }

    /// Return all conductor layers (Signal or Ground).
    pub fn conductor_layers(&self) -> Vec<&Layer> {
        self.layers
            .iter()
            .filter(|l| matches!(l.layer_type, LayerType::Signal | LayerType::Ground))
            .collect()
    }

    /// Return all dielectric layers.
    pub fn dielectric_layers(&self) -> Vec<&Layer> {
        self.layers
            .iter()
            .filter(|l| l.layer_type == LayerType::Dielectric)
            .collect()
    }

    /// Number of conductor layers.
    pub fn n_conductor_layers(&self) -> usize {
        self.layers
            .iter()
            .filter(|l| matches!(l.layer_type, LayerType::Signal | LayerType::Ground))
            .count()
    }

    /// Find a layer by name.
    pub fn find(&self, name: &str) -> Option<&Layer> {
        self.layers.iter().find(|l| l.name == name)
    }
}

impl Default for LayerStack {
    fn default() -> Self {
        Self::new("Untitled", "mm")
    }
}

impl Layer {
    pub fn new(name: impl Into<String>, layer_type: LayerType, thickness: f64) -> Self {
        Self {
            name: name.into(),
            layer_type,
            thickness,
            z_offset: 0.0,
            material: None,
            eps_r: None,
            loss_tangent: None,
            conductivity: None,
        }
    }

    /// Effective permittivity, looking up material if assigned (placeholder).
    pub fn effective_eps_r(&self, _materials: &[super::material::MaterialDef]) -> f64 {
        self.eps_r.unwrap_or(1.0)
    }

    /// Effective conductivity.
    pub fn effective_conductivity(&self) -> f64 {
        self.conductivity.unwrap_or(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_stack_zero_thickness() {
        let stack = LayerStack::new("test", "mm");
        assert!((stack.total_thickness() - 0.0).abs() < 1e-14);
    }

    #[test]
    fn four_layer_pcb() {
        let mut stack = LayerStack::new("4-layer PCB", "mm");
        stack.push_layer(Layer::new("Bottom Air", LayerType::Background, 10.0));
        stack.push_layer(Layer::new("Bottom Metal", LayerType::Signal, 0.035));
        stack.push_layer(Layer::new("Core", LayerType::Dielectric, 0.2));
        stack.push_layer(Layer::new("Inner 1", LayerType::Ground, 0.035));
        stack.push_layer(Layer::new("Prepreg", LayerType::Dielectric, 0.2));
        stack.push_layer(Layer::new("Top Metal", LayerType::Signal, 0.035));
        stack.push_layer(Layer::new("Top Air", LayerType::Background, 10.0));

        assert_eq!(stack.layers.len(), 7);
        assert_eq!(stack.n_conductor_layers(), 3);
        assert_eq!(stack.dielectric_layers().len(), 2);
        assert!(stack.total_thickness() > 20.0);
    }

    #[test]
    fn z_offset_auto_computed() {
        let mut stack = LayerStack::new("test", "mm");
        stack.push_layer(Layer::new("L1", LayerType::Dielectric, 1.0));
        stack.push_layer(Layer::new("L2", LayerType::Signal, 0.5));
        stack.push_layer(Layer::new("L3", LayerType::Dielectric, 1.0));

        assert!((stack.layers[0].z_offset - 0.0).abs() < 1e-14);
        assert!((stack.layers[1].z_offset - 1.0).abs() < 1e-14);
        assert!((stack.layers[2].z_offset - 1.5).abs() < 1e-14);
        assert!((stack.total_thickness() - 2.5).abs() < 1e-14);
    }
}

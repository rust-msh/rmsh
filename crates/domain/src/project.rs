use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::design::Design;
use crate::variable::{DatasetDefinition, Variable};

/// Top-level EMStudio project (the `.emsp` file).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmProject {
    pub metadata: ProjectMetadata,
    /// Project-level variables (`$` prefix in expressions).
    #[serde(default)]
    pub variables: HashMap<String, Variable>,
    /// Project-level datasets (frequency/temperature lookup tables).
    #[serde(default)]
    pub datasets: HashMap<String, DatasetDefinition>,
    #[serde(default)]
    pub designs: Vec<Design>,
}

impl Default for EmProject {
    fn default() -> Self {
        Self {
            metadata: ProjectMetadata::default(),
            variables: HashMap::new(),
            datasets: HashMap::new(),
            designs: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectMetadata {
    #[serde(default = "default_version")]
    pub version: String,
    #[serde(default = "default_application")]
    pub application: String,
    #[serde(default)]
    pub created_at: String,
    #[serde(default)]
    pub modified_at: String,
    #[serde(default)]
    pub author: String,
    #[serde(default)]
    pub description: String,
}

fn default_version() -> String {
    "1.0.0".to_string()
}
fn default_application() -> String {
    "EMStudio".to_string()
}

impl Default for ProjectMetadata {
    fn default() -> Self {
        Self {
            version: default_version(),
            application: default_application(),
            created_at: String::new(),
            modified_at: String::new(),
            author: String::new(),
            description: String::new(),
        }
    }
}

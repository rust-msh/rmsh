//! File I/O for .emsp project files (JSON format).
//!
//! Also provides lock file management and auto-save/recovery.

#[cfg(not(target_arch = "wasm32"))]
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::project::EmProject;
use crate::validation;

#[derive(Debug, Error)]
pub enum FileError {
    #[error("IO error: {0}")]
    Io(String),
    #[error("JSON parse error: {0}")]
    JsonParse(String),
    #[error("JSON write error: {0}")]
    JsonWrite(String),
    #[error("validation errors: {0:?}")]
    Validation(Vec<String>),
    #[error("project locked: {0}")]
    Locked(String),
    #[error("unsupported version: {0}")]
    UnsupportedVersion(String),
}

impl EmProject {
    /// Load a project from a JSON `.emsp` file.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load_json(path: &Path) -> Result<Self, FileError> {
        let content =
            std::fs::read_to_string(path).map_err(|e| FileError::Io(e.to_string()))?;
        Self::from_json_str(&content)
    }

    /// Parse a project from a JSON string.
    pub fn from_json_str(json: &str) -> Result<Self, FileError> {
        // Try to detect version for migration
        let raw: serde_json::Value =
            serde_json::from_str(json).map_err(|e| FileError::JsonParse(e.to_string()))?;

        let version = raw
            .get("metadata")
            .and_then(|m| m.get("version"))
            .and_then(|v| v.as_str())
            .unwrap_or("1.0.0");

        if !version.starts_with("1.") {
            return Err(FileError::UnsupportedVersion(version.to_string()));
        }

        let project: EmProject =
            serde_json::from_value(raw).map_err(|e| FileError::JsonParse(e.to_string()))?;

        Ok(project)
    }

    /// Save the project to a JSON `.emsp` file (pretty-printed).
    #[cfg(not(target_arch = "wasm32"))]
    pub fn save_json(&self, path: &Path) -> Result<(), FileError> {
        let json = self.to_json_string()?;
        std::fs::write(path, json).map_err(|e| FileError::Io(e.to_string()))
    }

    /// Serialize the project to a pretty-printed JSON string.
    pub fn to_json_string(&self) -> Result<String, FileError> {
        serde_json::to_string_pretty(self).map_err(|e| FileError::JsonWrite(e.to_string()))
    }

    /// Validate and return errors (empty = valid).
    pub fn validate(&self) -> Vec<String> {
        validation::validate_project(self)
            .into_iter()
            .map(|e| e.to_string())
            .collect()
    }

    /// Save to `.emsp.auto` for crash recovery.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn auto_save(&self, project_path: &Path) -> Result<(), FileError> {
        let auto_path = auto_save_path(project_path);
        self.save_json(&auto_path)
    }

    /// Recover from `.emsp.auto` if it exists.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn recover(project_path: &Path) -> Option<Self> {
        let auto_path = auto_save_path(project_path);
        if auto_path.exists() {
            Self::load_json(&auto_path).ok()
        } else {
            None
        }
    }

    /// Delete the auto-save file after a successful manual save.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn clear_auto_save(project_path: &Path) {
        let auto_path = auto_save_path(project_path);
        let _ = std::fs::remove_file(auto_path);
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn auto_save_path(project_path: &Path) -> PathBuf {
    project_path.with_extension("emsp.auto")
}

// ---------------------------------------------------------------------------
// Project lock (.emsp.lock)
// ---------------------------------------------------------------------------

/// File-based lock to prevent concurrent editing of a project.
#[cfg(not(target_arch = "wasm32"))]
pub struct ProjectLock {
    lock_path: PathBuf,
}

#[cfg(not(target_arch = "wasm32"))]
impl ProjectLock {
    /// Attempt to acquire an exclusive lock on the project file.
    pub fn acquire(project_path: &Path) -> Result<Self, FileError> {
        let lock_path = project_path.with_extension("emsp.lock");
        if lock_path.exists() {
            let content = std::fs::read_to_string(&lock_path).unwrap_or_default();
            return Err(FileError::Locked(format!(
                "Project is locked: {content}"
            )));
        }

        let info = serde_json::json!({
            "pid": std::process::id(),
            "locked_at": chrono_now_stub(),
        });

        std::fs::write(&lock_path, info.to_string())
            .map_err(|e| FileError::Io(e.to_string()))?;

        Ok(Self { lock_path })
    }

    /// Release the lock.
    pub fn release(&self) {
        let _ = std::fs::remove_file(&self.lock_path);
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl Drop for ProjectLock {
    fn drop(&mut self) {
        self.release();
    }
}

/// Stub for getting current time as a string (avoids chrono dependency).
fn chrono_now_stub() -> String {
    "2026-01-01T00:00:00Z".to_string()
}

// ---------------------------------------------------------------------------
// Result directory management
// ---------------------------------------------------------------------------

/// Create the result directory tree for a design.
///
/// ```text
/// {project_path}.results/
///   {design_id}/
///     {setup_name}/
///       mesh/
///       solutions/
///       fields/
///       ...
/// ```
#[cfg(not(target_arch = "wasm32"))]
pub fn create_result_dirs(
    project_path: &Path,
    design_id: &str,
    setup_names: &[&str],
) -> Result<PathBuf, FileError> {
    let results_dir = project_path.with_extension("emsp.results");
    let design_dir = results_dir.join(design_id);

    for setup_name in setup_names {
        let setup_dir = design_dir.join(setup_name);
        for subdir in &["mesh", "solutions", "fields"] {
            std::fs::create_dir_all(setup_dir.join(subdir))
                .map_err(|e| FileError::Io(e.to_string()))?;
        }
    }

    // Create exports dir
    std::fs::create_dir_all(design_dir.join("exports"))
        .map_err(|e| FileError::Io(e.to_string()))?;

    Ok(results_dir)
}

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Command: main thread -> worker
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WorkerCommand {
    /// Save a project to OPFS. `data` is MessagePack-encoded `Project`.
    SaveProject { data: Vec<u8> },
    /// Load a project by id from OPFS.
    LoadProject { id: String },
    /// Delete a project from OPFS.
    DeleteProject { id: String },
    /// List all projects stored in OPFS.
    ListProjects,
    /// Run the solver. `data` is MessagePack-encoded `Project`.
    Solve { data: Vec<u8> },
}

// ---------------------------------------------------------------------------
// Response: worker -> main thread
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WorkerResponse {
    /// Project saved successfully.
    SaveOk,
    /// Project data loaded. `data` is MessagePack-encoded `Project`.
    ProjectLoaded { data: Vec<u8> },
    /// Requested project was not found.
    ProjectNotFound { id: String },
    /// List of stored projects.
    ProjectList { entries: Vec<ProjectEntry> },
    /// Solve completed. `data` is MessagePack-encoded `SolveResult`.
    SolveResult { data: Vec<u8> },
    /// An error occurred in the worker.
    Error { message: String },
}

// ---------------------------------------------------------------------------
// Project index entry (lightweight metadata)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectEntry {
    pub id: String,
    pub title: String,
    /// Last modified timestamp in milliseconds since Unix epoch.
    pub modified_at: u64,
}

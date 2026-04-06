//! WASM backend that bridges the main thread to a Web Worker for OPFS + solver.
//!
//! This module is only compiled on `wasm32` targets.

use std::collections::HashMap;
use std::sync::mpsc;

use js_sys::{Array, Uint8Array};
use wasm_bindgen::prelude::*;

use emstudio_domain::worker_protocol::*;
use emstudio_domain::{Project, SolveResult};

use crate::{Backend, BackendError, RunMode};

// ---------------------------------------------------------------------------
// Worker status
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WorkerStatus {
    Starting,
    Ready,
    Error(String),
}

// ---------------------------------------------------------------------------
// WasmBackend
// ---------------------------------------------------------------------------

pub struct WasmBackend {
    worker: web_sys::Worker,
    status: WorkerStatus,
    response_rx: mpsc::Receiver<WorkerResponse>,
    // Keep the closures alive
    _onmessage: Closure<dyn FnMut(web_sys::MessageEvent)>,
    _onerror: Closure<dyn FnMut(web_sys::ErrorEvent)>,
    // Cached state
    projects_cache: HashMap<String, Project>,
    pending_loaded_project: Option<Project>,
    pending_solve_result: Option<SolveResult>,
    project_list: Vec<ProjectEntry>,
    // Error buffer for reporting to the UI
    last_error: Option<String>,
}

impl WasmBackend {
    /// Create a new WasmBackend and spawn the Web Worker.
    ///
    /// `worker_script_url` is the URL of the worker's JS bootstrap file,
    /// e.g. `"./worker_bootstrap.js"`.
    pub fn new(worker_script_url: &str) -> Result<Self, BackendError> {
        let worker = web_sys::Worker::new(worker_script_url)
            .map_err(|e| BackendError::IoError(format!("Failed to create Worker: {e:?}")))?;

        let (tx, rx) = mpsc::channel::<WorkerResponse>();

        // onmessage: deserialize WorkerResponse and send through channel
        let tx_msg = tx.clone();
        let onmessage = Closure::wrap(Box::new(move |event: web_sys::MessageEvent| {
            let data = event.data();
            let uint8 = Uint8Array::new(&data);
            let bytes = uint8.to_vec();
            match rmp_serde::from_slice::<WorkerResponse>(&bytes) {
                Ok(resp) => {
                    let _ = tx_msg.send(resp);
                }
                Err(e) => {
                    let _ = tx_msg.send(WorkerResponse::Error {
                        message: format!("Failed to deserialize worker response: {e}"),
                    });
                }
            }
        }) as Box<dyn FnMut(web_sys::MessageEvent)>);

        worker.set_onmessage(Some(onmessage.as_ref().unchecked_ref()));

        // onerror: report worker errors
        let tx_err = tx;
        let onerror = Closure::wrap(Box::new(move |event: web_sys::ErrorEvent| {
            let msg = event.message();
            let _ = tx_err.send(WorkerResponse::Error {
                message: format!("Worker error: {msg}"),
            });
        }) as Box<dyn FnMut(web_sys::ErrorEvent)>);

        worker.set_onerror(Some(onerror.as_ref().unchecked_ref()));

        Ok(Self {
            worker,
            status: WorkerStatus::Starting,
            response_rx: rx,
            _onmessage: onmessage,
            _onerror: onerror,
            projects_cache: HashMap::new(),
            pending_loaded_project: None,
            pending_solve_result: None,
            project_list: Vec::new(),
            last_error: None,
        })
    }

    /// Send a command to the Web Worker.
    fn send_command(&self, cmd: &WorkerCommand) -> Result<(), BackendError> {
        let data = rmp_serde::to_vec(cmd)
            .map_err(|e| BackendError::SerializeError(e.to_string()))?;
        let uint8 = Uint8Array::from(data.as_slice());
        let transfer = Array::of1(&uint8.buffer());
        self.worker
            .post_message_with_transfer(&uint8, &transfer)
            .map_err(|e| BackendError::IoError(format!("postMessage failed: {e:?}")))?;
        Ok(())
    }

    /// Get the current worker status.
    pub fn status(&self) -> &WorkerStatus {
        &self.status
    }

    /// Take the last error message, if any.
    pub fn take_last_error(&mut self) -> Option<String> {
        self.last_error.take()
    }

    /// Get the cached project list.
    pub fn project_list(&self) -> &[ProjectEntry] {
        &self.project_list
    }

    /// Request the worker to list all stored projects.
    pub fn request_project_list(&self) -> Result<(), BackendError> {
        self.send_command(&WorkerCommand::ListProjects)
    }
}

impl Backend for WasmBackend {
    fn save_project(&mut self, project: Project) -> Result<(), BackendError> {
        let data = rmp_serde::to_vec(&project)
            .map_err(|e| BackendError::SerializeError(e.to_string()))?;
        // Cache locally
        self.projects_cache.insert(project.id.clone(), project);
        self.send_command(&WorkerCommand::SaveProject { data })
    }

    fn load_project(&self, id: &str) -> Result<Project, BackendError> {
        // Return from cache if available
        if let Some(project) = self.projects_cache.get(id) {
            return Ok(project.clone());
        }
        // Otherwise, request from worker (result will arrive via poll)
        self.send_command(&WorkerCommand::LoadProject {
            id: id.to_string(),
        })?;
        Err(BackendError::ProjectNotFound(id.to_string()))
    }

    fn solve(&self, project: &Project) -> Result<SolveResult, BackendError> {
        let data = rmp_serde::to_vec(project)
            .map_err(|e| BackendError::SerializeError(e.to_string()))?;
        self.send_command(&WorkerCommand::Solve { data })?;
        // Result will arrive via poll
        Err(BackendError::IoError("solve pending".to_string()))
    }

    fn mode(&self) -> RunMode {
        RunMode::LocalFirst
    }

    fn poll(&mut self) {
        while let Ok(resp) = self.response_rx.try_recv() {
            match resp {
                WorkerResponse::SaveOk => {
                    if self.status == WorkerStatus::Starting {
                        self.status = WorkerStatus::Ready;
                    }
                }
                WorkerResponse::ProjectLoaded { data } => {
                    match rmp_serde::from_slice::<Project>(&data) {
                        Ok(project) => {
                            self.pending_loaded_project = Some(project.clone());
                            self.projects_cache.insert(project.id.clone(), project);
                        }
                        Err(e) => {
                            self.last_error =
                                Some(format!("Failed to deserialize loaded project: {e}"));
                        }
                    }
                }
                WorkerResponse::ProjectNotFound { id } => {
                    self.last_error = Some(format!("Project not found: {id}"));
                }
                WorkerResponse::ProjectList { entries } => {
                    self.project_list = entries;
                }
                WorkerResponse::SolveResult { data } => {
                    match rmp_serde::from_slice::<SolveResult>(&data) {
                        Ok(result) => {
                            self.pending_solve_result = Some(result);
                        }
                        Err(e) => {
                            self.last_error =
                                Some(format!("Failed to deserialize solve result: {e}"));
                        }
                    }
                }
                WorkerResponse::Error { message } => {
                    self.status = WorkerStatus::Error(message.clone());
                    self.last_error = Some(message);
                }
            }
        }
    }

    fn take_solve_result(&mut self) -> Option<SolveResult> {
        self.pending_solve_result.take()
    }

    fn take_loaded_project(&mut self) -> Option<Project> {
        self.pending_loaded_project.take()
    }
}

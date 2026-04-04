//! Web Worker entry point for EMStudio's Local-First mode.
//!
//! This crate compiles to a separate WASM module that runs inside a Web Worker.
//! It owns the OPFS file system handles and the solver, communicating with the
//! main thread via `postMessage` using MessagePack-encoded commands/responses.

mod opfs;

use js_sys::{Array, Uint8Array};
use wasm_bindgen::prelude::*;
use wasm_bindgen_futures::spawn_local;

use emstudio_domain::worker_protocol::*;
use emstudio_domain::Project;
use emstudio_solver::{PlaceholderSolver, Solver};

// ---------------------------------------------------------------------------
// Helper: get DedicatedWorkerGlobalScope
// ---------------------------------------------------------------------------

fn worker_global_scope() -> web_sys::DedicatedWorkerGlobalScope {
    js_sys::global()
        .dyn_into::<web_sys::DedicatedWorkerGlobalScope>()
        .expect("not running in a DedicatedWorkerGlobalScope")
}

// ---------------------------------------------------------------------------
// Helper: post a WorkerResponse back to the main thread
// ---------------------------------------------------------------------------

fn post_response(resp: &WorkerResponse) {
    let data = rmp_serde::to_vec(resp).expect("failed to serialize WorkerResponse");
    let uint8 = Uint8Array::from(data.as_slice());
    let scope = worker_global_scope();
    // Transfer the underlying ArrayBuffer for zero-copy
    let transfer = Array::of1(&uint8.buffer());
    scope
        .post_message_with_transfer(&uint8, &transfer)
        .expect("post_message failed");
}

// ---------------------------------------------------------------------------
// Command dispatch
// ---------------------------------------------------------------------------

async fn handle_command(cmd: WorkerCommand) {
    let resp = match cmd {
        WorkerCommand::SaveProject { data } => handle_save_project(data).await,
        WorkerCommand::LoadProject { id } => handle_load_project(id).await,
        WorkerCommand::DeleteProject { id } => handle_delete_project(id).await,
        WorkerCommand::ListProjects => handle_list_projects().await,
        WorkerCommand::Solve { data } => handle_solve(data),
    };
    post_response(&resp);
}

async fn handle_save_project(data: Vec<u8>) -> WorkerResponse {
    // Deserialize just to get the project id and title for the index
    let project: Project = match rmp_serde::from_slice(&data) {
        Ok(p) => p,
        Err(e) => {
            return WorkerResponse::Error {
                message: format!("Failed to deserialize project for save: {e}"),
            }
        }
    };

    let dir = match opfs::ensure_projects_dir().await {
        Ok(d) => d,
        Err(e) => {
            return WorkerResponse::Error {
                message: format!("OPFS ensure_projects_dir failed: {e:?}"),
            }
        }
    };

    let filename = format!("{}.emsp", project.id);
    if let Err(e) = opfs::write_file(&dir, &filename, &data).await {
        return WorkerResponse::Error {
            message: format!("OPFS write failed: {e:?}"),
        };
    }

    // Update index
    if let Err(e) = update_index(&dir, &project).await {
        // Non-fatal: project is saved, index update failed
        web_sys::console::warn_1(&format!("Index update failed: {e:?}").into());
    }

    WorkerResponse::SaveOk
}

async fn handle_load_project(id: String) -> WorkerResponse {
    let dir = match opfs::ensure_projects_dir().await {
        Ok(d) => d,
        Err(e) => {
            return WorkerResponse::Error {
                message: format!("OPFS ensure_projects_dir failed: {e:?}"),
            }
        }
    };

    let filename = format!("{id}.emsp");
    match opfs::read_file(&dir, &filename).await {
        Ok(data) => WorkerResponse::ProjectLoaded { data },
        Err(_) => WorkerResponse::ProjectNotFound { id },
    }
}

async fn handle_delete_project(id: String) -> WorkerResponse {
    let dir = match opfs::ensure_projects_dir().await {
        Ok(d) => d,
        Err(e) => {
            return WorkerResponse::Error {
                message: format!("OPFS ensure_projects_dir failed: {e:?}"),
            }
        }
    };

    let filename = format!("{id}.emsp");
    if let Err(e) = opfs::delete_file(&dir, &filename).await {
        return WorkerResponse::Error {
            message: format!("OPFS delete failed: {e:?}"),
        };
    }

    // Remove from index
    if let Err(e) = remove_from_index(&dir, &id).await {
        web_sys::console::warn_1(&format!("Index removal failed: {e:?}").into());
    }

    WorkerResponse::SaveOk
}

async fn handle_list_projects() -> WorkerResponse {
    let dir = match opfs::ensure_projects_dir().await {
        Ok(d) => d,
        Err(e) => {
            return WorkerResponse::Error {
                message: format!("OPFS ensure_projects_dir failed: {e:?}"),
            }
        }
    };

    match load_index(&dir).await {
        Ok(entries) => WorkerResponse::ProjectList { entries },
        Err(e) => WorkerResponse::Error {
            message: format!("Failed to load project index: {e:?}"),
        },
    }
}

fn handle_solve(data: Vec<u8>) -> WorkerResponse {
    let project: Project = match rmp_serde::from_slice(&data) {
        Ok(p) => p,
        Err(e) => {
            return WorkerResponse::Error {
                message: format!("Failed to deserialize project for solve: {e}"),
            }
        }
    };

    let solver = PlaceholderSolver;
    let result = solver.solve(&project.model);

    match rmp_serde::to_vec(&result) {
        Ok(result_data) => WorkerResponse::SolveResult { data: result_data },
        Err(e) => WorkerResponse::Error {
            message: format!("Failed to serialize solve result: {e}"),
        },
    }
}

// ---------------------------------------------------------------------------
// Project index management (index.json in OPFS)
// ---------------------------------------------------------------------------

const INDEX_FILE: &str = "index.json";

async fn load_index(dir: &JsValue) -> Result<Vec<ProjectEntry>, String> {
    match opfs::read_file(dir, INDEX_FILE).await {
        Ok(data) => serde_json::from_slice(&data).map_err(|e| format!("parse index: {e}")),
        Err(_) => Ok(Vec::new()), // No index file yet — return empty list
    }
}

async fn save_index(dir: &JsValue, entries: &[ProjectEntry]) -> Result<(), String> {
    let data = serde_json::to_vec(entries).map_err(|e| format!("serialize index: {e}"))?;
    opfs::write_file(dir, INDEX_FILE, &data)
        .await
        .map_err(|e| format!("write index: {e:?}"))
}

async fn update_index(dir: &JsValue, project: &Project) -> Result<(), String> {
    let mut entries = load_index(dir).await?;

    // Get current time in millis
    let now = js_sys::Date::now() as u64;

    // Update or insert entry
    if let Some(entry) = entries.iter_mut().find(|e| e.id == project.id) {
        entry.title = project.title.clone();
        entry.modified_at = now;
    } else {
        entries.push(ProjectEntry {
            id: project.id.clone(),
            title: project.title.clone(),
            modified_at: now,
        });
    }

    save_index(dir, &entries).await
}

async fn remove_from_index(dir: &JsValue, id: &str) -> Result<(), String> {
    let mut entries = load_index(dir).await?;
    entries.retain(|e| e.id != id);
    save_index(dir, &entries).await
}

// ---------------------------------------------------------------------------
// Worker entry point
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn worker_entry() {
    console_error_panic_hook::set_once();

    let scope = worker_global_scope();

    // Register the onmessage handler
    let onmessage = Closure::wrap(Box::new(move |event: web_sys::MessageEvent| {
        let data = event.data();
        let uint8 = Uint8Array::new(&data);
        let bytes = uint8.to_vec();

        match rmp_serde::from_slice::<WorkerCommand>(&bytes) {
            Ok(cmd) => {
                spawn_local(async move {
                    handle_command(cmd).await;
                });
            }
            Err(e) => {
                let resp = WorkerResponse::Error {
                    message: format!("Failed to deserialize command: {e}"),
                };
                post_response(&resp);
            }
        }
    }) as Box<dyn FnMut(web_sys::MessageEvent)>);

    scope.set_onmessage(Some(onmessage.as_ref().unchecked_ref()));
    onmessage.forget(); // Leak the closure to keep it alive

    // Signal that the worker is ready
    let ready = WorkerResponse::SaveOk; // Reuse as a "ready" signal
    post_response(&ready);
}

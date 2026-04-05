# EmStudio Storage & File System Backends Analysis

**Updated**: April 4, 2026 | **Status**: Comprehensive exploration complete

---

## Executive Summary

### Current State
- **Desktop**: ✅ File-based storage via MessagePack (.emsp files)
- **Web**: ❌ NO persistent storage (everything lost on page reload)
- **Worker usage**: ❌ NOT IMPLEMENTED
- **OPFS support**: ❌ NOT IMPLEMENTED
- **Auto-save**: ❌ NOT IMPLEMENTED
- **Version control**: ❌ NOT IMPLEMENTED

### Key Finding
The codebase has **clean separation** between:
- **File I/O layer** (`crates/infra/src/lib.rs`): `save_project_to_file()`, `load_project_from_file()`
- **Backend abstraction** (`Backend` trait): Storage-agnostic interface
- **App layer** (`crates/app/src/lib.rs`): UI-driven file operations

**BUT**: File I/O is **NOT integrated** into the Backend trait - this is a critical gap.

---

## 1. Current File I/O Implementation

### Location
`/crates/infra/src/lib.rs` (lines 32-48)

### Code
```rust
pub fn save_project_to_file(project: &Project, path: &Path) -> Result<(), BackendError> {
    let data = rmp_serde::to_vec(project)
        .map_err(|e| BackendError::SerializeError(e.to_string()))?;
    std::fs::write(path, data)
        .map_err(|e| BackendError::IoError(e.to_string()))?;
    Ok(())
}

pub fn load_project_from_file(path: &Path) -> Result<Project, BackendError> {
    let data = std::fs::read(path)
        .map_err(|e| BackendError::IoError(e.to_string()))?;
    let project: Project = rmp_serde::from_slice(&data)
        .map_err(|e| BackendError::DeserializeError(e.to_string()))?;
    Ok(project)
}
```

### File Format
- **Name**: `.emsp` (EmStudio Project)
- **Encoding**: MessagePack (binary)
- **Structure**: Serialized `Project` struct
- **Serialization**: serde (auto-derived via `#[derive(Serialize, Deserialize)]`)

### Advantages
- ✅ **Fast**: Binary format, no parsing overhead
- ✅ **Compact**: Smaller files than JSON
- ✅ **Type-safe**: Rust serde ensures correctness
- ✅ **Cross-platform**: Works on Windows/macOS/Linux
- ✅ **WASM-compatible**: MessagePack works in browsers (with proper file access)

### Limitations
- ❌ **Not human-readable**: Can't inspect with text editors
- ❌ **No versioning**: Breaking changes require migration code
- ❌ **Not design-compliant**: Spec calls for JSON in `em-project-file-design.md`
- ❌ **No compression**: Large projects = large files
- ❌ **No incremental save**: Full file rewrite every time
- ❌ **Not integrated in Backend trait**: Loose coupling

---

## 2. Backend Abstraction Layer

### Location
`/crates/infra/src/lib.rs` (lines 54-146)

### Backend Trait
```rust
pub trait Backend {
    fn save_project(&mut self, project: Project) -> Result<(), BackendError>;
    fn load_project(&self, id: &str) -> Result<Project, BackendError>;
    fn solve(&self, project: &Project) -> Result<SolveResult, BackendError>;
    fn mode(&self) -> RunMode;
}
```

### Implementations

#### StandaloneBackend (In-Memory)
```rust
pub struct StandaloneBackend {
    projects: HashMap<String, Project>,
    solver: PlaceholderSolver,
}

impl Backend for StandaloneBackend {
    fn save_project(&mut self, project: Project) -> Result<(), BackendError> {
        self.projects.insert(project.id.clone(), project);
        Ok(())
    }
    fn load_project(&self, id: &str) -> Result<Project, BackendError> {
        self.projects.get(id).cloned()
            .ok_or_else(|| BackendError::ProjectNotFound(id.to_string()))
    }
    // ...
}
```

**Characteristics**:
- ✅ No persistence (data lost on exit)
- ✅ Projects stored in RAM (fast access)
- ✅ Multi-project support (HashMap)
- ❌ No file system interaction
- ❌ No resource cleanup

#### CloudBackend (Remote API - STUB)
```rust
pub struct CloudBackend {
    endpoint: String,
}

impl Backend for CloudBackend {
    fn save_project(&mut self, _project: Project) -> Result<(), BackendError> {
        Ok(())  // ← No-op!
    }
    // ...
}
```

**Characteristics**:
- ❌ NOT IMPLEMENTED (all methods are stubs)
- ❌ No actual HTTP communication
- ❌ No error handling

### Critical Gap
**File I/O functions are NOT methods on Backend trait**:
```rust
// Current: App calls these directly
save_project_to_file(&self.project, path)?;    // From infra module
load_project_from_file(path)?;

// Should be: App calls through Backend trait
backend.save_project(project)?;
backend.load_project(id)?;
```

This means:
- Backend abstraction is **incomplete**
- File operations are **not pluggable** per backend
- **Hard to test** different storage implementations
- **Prevents OPFS** integration (would need Backend method)

---

## 3. App Layer Integration

### Location
`/crates/app/src/lib.rs`

### File Operations Methods
```rust
pub fn save_to(&mut self, path: &std::path::Path) {
    match save_project_to_file(&self.project, path) {
        Ok(()) => {
            self.current_file = Some(path.to_path_buf());
            self.unsaved_changes = false;
            self.status_text = format!("Saved: {}", path.display());
        }
        Err(e) => { /* error handling */ }
    }
}

pub fn open_from(&mut self, path: &std::path::Path) {
    match load_project_from_file(path) {
        Ok(project) => {
            self.project = project;
            self.current_file = Some(path.to_path_buf());
            self.unsaved_changes = false;
        }
        Err(e) => { /* error handling */ }
    }
}
```

### Async File Dialog
```rust
fn spawn_open_dialog(&self, ctx: &egui::Context) {
    let tx = self.file_dialog_tx.clone();
    spawn_future(async move {
        let file = rfd::AsyncFileDialog::new()
            .add_filter("EmStudio Project", &["emsp"])
            .pick_file()
            .await;
        if let Some(handle) = file {
            #[cfg(not(target_arch = "wasm32"))]
            let path = handle.path().to_path_buf();
            #[cfg(target_arch = "wasm32")]
            let path = PathBuf::from(handle.file_name());

            let _ = tx.send(FileDialogResult::OpenFile(path));
            ctx.request_repaint();
        }
    });
}
```

**Characteristics**:
- ✅ Cross-platform (native + WASM dialogs)
- ✅ Async (non-blocking)
- ⚠️ WASM returns filename only (no actual file access)
- ❌ No progress feedback
- ❌ No resume on partial load

### State Management
```rust
pub struct App {
    current_file: Option<PathBuf>,
    unsaved_changes: bool,
    status_text: String,
    // ...
}
```

**Limitations**:
- ❌ No auto-save timer
- ❌ No backup/recovery mechanism
- ❌ No file locking (for multi-user scenarios)
- ❌ No undo/redo system

---

## 4. Missing from Design Spec

### File Formats (from docs/)
Design calls for TWO formats:
1. **JSON** (human-readable, versioning support)
2. **Binary** (optional, for performance)

Current implementation:
- ❌ JSON not used
- ✅ Binary (MessagePack) only

### Auto-Save System
**Design**: `project.emsp` + `project.emsp.auto` (backup)

**Current**: ❌ NOT IMPLEMENTED

### File Locking
**Design**: `project.emsp.lock` (prevent concurrent writes)

**Current**: ❌ NOT IMPLEMENTED

### Result Directory Structure
**Design**: 
```
Project.emsp.results/
├── design-001/
│   ├── validation_report.json
│   ├── Setup1/
│   │   ├── convergence.json
│   │   ├── mesh_stats.json
│   │   ├── solutions/
│   │   └── fields/
```

**Current**: ❌ NOT IMPLEMENTED

---

## 5. OPFS (Origin Private File System) - Web Storage

### Design (from docs/em-feature-design-and-progress.md, §1.5.1)

```
┌─────────────────────────────────────────────────────┐
│                 Browser Main Thread                  │
│         egui (WASM) + WebGPU Rendering              │
├─────────────────────────────────────────────────────┤
│              Web Worker (Backend Thread)             │
│   ┌─────────────────┐  ┌──────────────────────┐    │
│   │ Rem Solver      │  │ OPFS Storage         │    │
│   │ (WASM compiled) │  │ ├── project.emsp     │    │
│   │                 │  │ ├── results/         │    │
│   │                 │  │ └── materials.emsm   │    │
│   └─────────────────┘  └──────────────────────┘    │
└─────────────────────────────────────────────────────┘
```

### What is OPFS?
- **Origin Private File System** (web standard)
- Persistent storage in browser (like `localStorage` but file-based)
- NOT accessible by JavaScript (only via FileSystemHandle)
- Survives page reloads and browser restarts
- Separate sandbox per origin

### Current Status
- ❌ NOT IMPLEMENTED
- 🔲 Milestone 10 task (Platform & Deployment)
- No Web Worker setup
- No OPFS API integration

### Implementation Requirements
1. **Main thread**:
   - Spawn Web Worker
   - Send/receive messages for solver progress
   - Handle UI updates

2. **Web Worker**:
   - Access OPFS (via FileSystemHandle)
   - Run Rem solver (WASM)
   - Read/write project files
   - Stream results back to main thread

3. **OPFS API** (browser-side):
   ```javascript
   // Main thread
   const root = await navigator.storage.getDirectory();
   const projectFile = await root.getFileHandle('project.emsp', { create: true });
   
   // Worker thread
   const writable = await handle.createWritable();
   await writable.write(new Uint8Array([...]));
   await writable.close();
   ```

---

## 6. Web Worker Architecture (PLANNED)

### Current Async Pattern
```rust
fn spawn_future<F: std::future::Future<Output = ()> + 'static>(f: F) {
    #[cfg(target_arch = "wasm32")]
    wasm_bindgen_futures::spawn_local(f);  // ← Only spawn_local!
    #[cfg(not(target_arch = "wasm32"))]
    pollster::block_on(f);
}
```

**Problem**: `spawn_local()` runs on main thread, can block UI

### Proposed Web Worker Usage
```rust
// Main thread (app)
let worker = Worker::new("solver_worker.js")?;
worker.post_message(&serde_wasm_bindgen::to_value(&project)?);

// Worker thread (solver_worker.rs)
#[wasm_bindgen]
pub fn solve_in_worker(project_json: JsValue) {
    let project: Project = serde_wasm_bindgen::from_value(project_json)?;
    let result = solver.solve(&project);
    web_sys::window()
        .postMessage(&serde_wasm_bindgen::to_value(&result)?);
}

// Main thread (receive)
worker.set_onmessage(Some(Box::new(|msg| {
    let result: SolveResult = serde_wasm_bindgen::from_value(msg.data())?;
    app.project.last_result = Some(result);
})));
```

### Key Limitations
- ❌ No native `web_sys::Worker` in current codebase
- ❌ Requires separate `.js` worker file (or WASM blob)
- ❌ Message serialization overhead
- ⚠️ Worker doesn't share memory (must copy projects)

---

## 7. Serialization Details

### Current (MessagePack via serde)
```rust
#[derive(Serialize, Deserialize)]
pub struct Project {
    pub id: String,
    pub title: String,
    pub model: EmModel,
    pub status: SimulationStatus,
    pub last_result: Option<SolveResult>,
}

// Save
let bytes = rmp_serde::to_vec(&project)?;
std::fs::write("file.emsp", bytes)?;

// Load
let bytes = std::fs::read("file.emsp")?;
let project: Project = rmp_serde::from_slice(&bytes)?;
```

**Serde Stack**:
1. `serde` (trait definitions)
2. `rmp-serde` (MessagePack codec)
3. Custom derives (`#[derive(Serialize, Deserialize)]`)

### Design Spec (JSON)
```json
{
  "version": "1.0",
  "project": {
    "id": "proj-001",
    "title": "My Design",
    "designs": [
      {
        "id": "design-001",
        "name": "HFSS Analysis",
        "geometry": { ... },
        "materials": { ... },
        "boundaries": { ... }
      }
    ]
  }
}
```

**Mismatch**:
- ❌ Design calls for multi-design support (not in current Domain)
- ❌ Design calls for JSON (not MessagePack)
- ❌ Design calls for versioning (no migration strategy)

---

## 8. Cross-Platform File Access

### Desktop (Native)
```rust
// Full filesystem access
std::fs::write(path, bytes)?;  // ✅ Works everywhere
std::fs::read(path)?;

// Native dialog
rfd::AsyncFileDialog::new()
    .pick_file()
    .await  // ✅ Returns real FileHandle
```

### Web (WASM)
```rust
// NO filesystem access
std::fs::write(path, bytes)  // ❌ Not available in browser

// WASM dialog (limited)
rfd::AsyncFileDialog::new()
    .pick_file()
    .await  // ⚠️ Returns FileHandle, but only for user-selected files

// Browser API access (needs OPFS implementation)
#[cfg(target_arch = "wasm32")]
async fn opfs_write(path: &str, bytes: &[u8]) {
    // Requires Web Worker + FileSystemHandle
    // NOT IMPLEMENTED YET
}
```

### Conditional Compilation
```rust
#[cfg(not(target_arch = "wasm32"))]
fn code_for_desktop() { 
    std::fs::write(...)?;  // ✅ OK
}

#[cfg(target_arch = "wasm32")]
fn code_for_web() {
    // Use OPFS or IndexedDB
    // Currently: NOTHING (data lost on reload)
}
```

---

## 9. Testing Infrastructure

### Test Files
- `/crates/infra/tests/` (none currently)
- `/crates/app/tests/e2e.rs` (integration tests)
- Domain unit tests in `/crates/domain/src/lib.rs`

### Coverage
```rust
#[test]
fn save_and_load_roundtrip() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("test.emsp");

    let original = sample_project();
    save_project_to_file(&original, &path).unwrap();

    let loaded = load_project_from_file(&path).unwrap();
    assert_eq!(loaded.title, original.title);
}
```

**What's tested**:
- ✅ File I/O roundtrip
- ✅ Serialization correctness
- ✅ Error cases (corrupt files, missing files)

**What's NOT tested**:
- ❌ Concurrent file access
- ❌ Large file performance
- ❌ OPFS operations
- ❌ Web Worker communication
- ❌ Auto-save/recovery

---

## 10. Key Files Summary

| File | Purpose | Lines | Status |
|------|---------|-------|--------|
| `/crates/infra/src/lib.rs` | File I/O functions | 255 | ⭐⭐⭐ Complete |
| `/crates/app/src/lib.rs` | App layer (UI integration) | 545 | ⭐⭐⭐ Complete |
| `/crates/domain/src/lib.rs` | Domain models | 141 | ⭐⭐ Incomplete |
| `/crates/render/src/lib.rs` | 3D viewport | 150+ | ⭐⭐⭐⭐⭐ Complete |
| `/crates/main/src/main.rs` | Native entry point | 52 | ⭐⭐⭐ Complete |
| `/crates/main/src/lib.rs` | WASM entry point | 37 | ⭐⭐ Not tested |
| `/crates/app/tests/e2e.rs` | Integration tests | 220 | ⭐⭐⭐ Good |
| `docs/em-project-file-design.md` | File format spec | 3900 | ⭐⭐⭐⭐⭐ Detailed |

---

## 11. Recommendations

### Short-term (Implement next)
1. **Integrate file I/O into Backend trait**
   - Add `save_project()` and `load_project()` methods
   - Remove direct filesystem calls from App layer
   - Enable pluggable storage backends

2. **Add OPFS backend stub**
   - New `struct OpfsBackend`
   - Placeholder methods (return errors for now)
   - Framework for Milestone 10

3. **JSON format support**
   - Add optional JSON serialization alongside MessagePack
   - Version field for migration support
   - Human-readable for debugging

### Medium-term (Milestone 3-5)
1. **Auto-save system**
   - Timer-based periodic save
   - `.emsp.auto` backup file
   - Crash recovery on startup

2. **File locking**
   - `.emsp.lock` metadata file
   - Concurrent access prevention

3. **Result directory**
   - Structured results/ subdirectory
   - JSON metadata for each solve
   - Field data storage (.emsfld)

### Long-term (Milestone 10)
1. **Web Worker implementation**
   - Solver runs off main thread
   - OPFS file access from worker
   - Progress streaming to UI

2. **OPFS backend (full)**
   - File operations via FileSystemHandle
   - Persistent storage in browser
   - Offline-capable app

3. **Cloud backend (production)**
   - HTTP API integration
   - Multi-user session management
   - Server-side file storage

---

## 12. Technical Debt

| Issue | Impact | Fix Effort | Priority |
|-------|--------|-----------|----------|
| File I/O outside Backend trait | Blocks storage abstraction | Medium | HIGH |
| No JSON format support | Design compliance | Medium | MEDIUM |
| No domain model expansion | Blocks M3, M4, M5 | Large | HIGH |
| No WASM storage layer | Web deployment broken | Large | HIGH |
| No auto-save | Data loss risk | Small | MEDIUM |
| No file versioning | Migration required | Medium | LOW |
| No compression | Large file bloat | Small | LOW |

---

## Conclusion

EmStudio has a **functional desktop file system** but **no web storage** and **incomplete backend abstraction**. The critical path to web deployment requires:

1. **Refactor Backend trait** to include file I/O methods
2. **Implement OPFS backend** with Web Worker support
3. **Expand domain model** to match design spec

The **rendering engine is production-ready** but depends on the above for complete functionality.


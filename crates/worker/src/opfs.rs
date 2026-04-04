//! OPFS (Origin Private File System) wrapper functions.
//!
//! Uses `js_sys::Reflect` for OPFS APIs that lack stable `web-sys` bindings.

use js_sys::{Array, Object, Reflect, Uint8Array};
use wasm_bindgen::prelude::*;
use wasm_bindgen_futures::JsFuture;

// ---------------------------------------------------------------------------
// Directory handles
// ---------------------------------------------------------------------------

/// Get the OPFS root directory via `navigator.storage.getDirectory()`.
pub async fn get_root_dir() -> Result<JsValue, JsValue> {
    let global: JsValue = js_sys::global().into();
    let navigator = Reflect::get(&global, &"navigator".into())?;
    let storage = Reflect::get(&navigator, &"storage".into())?;
    let promise = Reflect::apply(
        &Reflect::get(&storage, &"getDirectory".into())?.into(),
        &storage,
        &Array::new(),
    )?;
    JsFuture::from(js_sys::Promise::from(promise)).await
}

/// Get or create a subdirectory inside the given parent directory handle.
pub async fn get_or_create_subdir(parent: &JsValue, name: &str) -> Result<JsValue, JsValue> {
    let opts = Object::new();
    Reflect::set(&opts, &"create".into(), &JsValue::TRUE)?;
    let args = Array::of2(&name.into(), &opts.into());
    let get_dir_fn = Reflect::get(parent, &"getDirectoryHandle".into())?;
    let promise = Reflect::apply(&get_dir_fn.into(), parent, &args)?;
    JsFuture::from(js_sys::Promise::from(promise)).await
}

/// Ensure the `/emstudio/projects/` directory chain exists in OPFS.
pub async fn ensure_projects_dir() -> Result<JsValue, JsValue> {
    let root = get_root_dir().await?;
    let emstudio = get_or_create_subdir(&root, "emstudio").await?;
    get_or_create_subdir(&emstudio, "projects").await
}

// ---------------------------------------------------------------------------
// File operations
// ---------------------------------------------------------------------------

/// Write `data` bytes to a file named `name` in the given directory handle.
pub async fn write_file(dir: &JsValue, name: &str, data: &[u8]) -> Result<(), JsValue> {
    // dir.getFileHandle(name, { create: true })
    let opts = Object::new();
    Reflect::set(&opts, &"create".into(), &JsValue::TRUE)?;
    let args = Array::of2(&name.into(), &opts.into());
    let get_fh_fn = Reflect::get(dir, &"getFileHandle".into())?;
    let file_handle =
        JsFuture::from(js_sys::Promise::from(Reflect::apply(&get_fh_fn.into(), dir, &args)?))
            .await?;

    // fileHandle.createWritable()
    let create_writable_fn = Reflect::get(&file_handle, &"createWritable".into())?;
    let writable = JsFuture::from(js_sys::Promise::from(Reflect::apply(
        &create_writable_fn.into(),
        &file_handle,
        &Array::new(),
    )?))
    .await?;

    // writable.write(data)
    let uint8 = Uint8Array::from(data);
    let write_fn = Reflect::get(&writable, &"write".into())?;
    let write_args = Array::of1(&uint8.into());
    JsFuture::from(js_sys::Promise::from(Reflect::apply(
        &write_fn.into(),
        &writable,
        &write_args,
    )?))
    .await?;

    // writable.close()
    let close_fn = Reflect::get(&writable, &"close".into())?;
    JsFuture::from(js_sys::Promise::from(Reflect::apply(
        &close_fn.into(),
        &writable,
        &Array::new(),
    )?))
    .await?;

    Ok(())
}

/// Read file contents as `Vec<u8>` from a file named `name` in the given directory.
pub async fn read_file(dir: &JsValue, name: &str) -> Result<Vec<u8>, JsValue> {
    // dir.getFileHandle(name)
    let args = Array::of1(&name.into());
    let get_fh_fn = Reflect::get(dir, &"getFileHandle".into())?;
    let file_handle =
        JsFuture::from(js_sys::Promise::from(Reflect::apply(&get_fh_fn.into(), dir, &args)?))
            .await?;

    // fileHandle.getFile()
    let get_file_fn = Reflect::get(&file_handle, &"getFile".into())?;
    let file = JsFuture::from(js_sys::Promise::from(Reflect::apply(
        &get_file_fn.into(),
        &file_handle,
        &Array::new(),
    )?))
    .await?;

    // file.arrayBuffer()
    let array_buf_fn = Reflect::get(&file, &"arrayBuffer".into())?;
    let buf = JsFuture::from(js_sys::Promise::from(Reflect::apply(
        &array_buf_fn.into(),
        &file,
        &Array::new(),
    )?))
    .await?;

    let uint8 = Uint8Array::new(&buf);
    Ok(uint8.to_vec())
}

/// Delete a file named `name` from the given directory handle.
pub async fn delete_file(dir: &JsValue, name: &str) -> Result<(), JsValue> {
    let args = Array::of1(&name.into());
    let remove_fn = Reflect::get(dir, &"removeEntry".into())?;
    JsFuture::from(js_sys::Promise::from(Reflect::apply(
        &remove_fn.into(),
        dir,
        &args,
    )?))
    .await?;
    Ok(())
}

/// List all file names in the given directory handle.
pub async fn list_files(dir: &JsValue) -> Result<Vec<String>, JsValue> {
    // Use dir.keys() async iterator
    let keys_fn = Reflect::get(dir, &"keys".into())?;
    let iterator = Reflect::apply(&keys_fn.into(), dir, &Array::new())?;

    let mut names = Vec::new();
    loop {
        let next_fn = Reflect::get(&iterator, &"next".into())?;
        let result = JsFuture::from(js_sys::Promise::from(Reflect::apply(
            &next_fn.into(),
            &iterator,
            &Array::new(),
        )?))
        .await?;

        let done = Reflect::get(&result, &"done".into())?;
        if done.as_bool().unwrap_or(true) {
            break;
        }

        let value = Reflect::get(&result, &"value".into())?;
        if let Some(s) = value.as_string() {
            names.push(s);
        }
    }

    Ok(names)
}

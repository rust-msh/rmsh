#[cfg(target_arch = "wasm32")]
mod wasm_app {
    use eframe::WebRunner;
    use wasm_bindgen::JsCast;
    use wasm_bindgen::prelude::*;

    use emstudio_app::App;
    use emstudio_domain::Edition;
    use emstudio_infra::RunMode;

    async fn start_on_canvas(canvas_id: &str) -> Result<(), JsValue> {
        console_error_panic_hook::set_once();

        let window =
            web_sys::window().ok_or_else(|| JsValue::from_str("window is unavailable"))?;
        let document = window
            .document()
            .ok_or_else(|| JsValue::from_str("document is unavailable"))?;
        let canvas = document
            .get_element_by_id(canvas_id)
            .ok_or_else(|| JsValue::from_str("missing canvas"))?
            .dyn_into::<web_sys::HtmlCanvasElement>()?;

        WebRunner::new()
            .start(
                canvas,
                eframe::WebOptions::default(),
                Box::new(|cc| Ok(Box::new(App::new(RunMode::LocalFirst, Edition::Basic, cc)))),
            )
            .await
    }

    #[wasm_bindgen]
    pub async fn start_with_canvas_id(canvas_id: &str) -> Result<(), JsValue> {
        start_on_canvas(canvas_id).await
    }

    #[wasm_bindgen]
    pub async fn start_with_canvas_and_project(
        canvas_id: &str,
        _project_id: &str,
        _project_name: &str,
        _project_description: &str,
        _owner_user_id: &str,
        _member_count: u32,
    ) -> Result<(), JsValue> {
        start_on_canvas(canvas_id).await
    }
}

#[cfg(not(target_arch = "wasm32"))]
pub fn host_build_marker() {
    // Intentionally empty: allows host cargo check for this crate.
}

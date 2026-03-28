#[cfg(target_arch = "wasm32")]
mod wasm_app {
    use eframe::WebRunner;
    use wasm_bindgen::JsCast;
    use wasm_bindgen::prelude::*;

    use emstudio_app::EmStudioApp;

    #[wasm_bindgen(start)]
    pub async fn start() -> Result<(), JsValue> {
        console_error_panic_hook::set_once();

        let window =
            web_sys::window().ok_or_else(|| JsValue::from_str("window is unavailable"))?;
        let document = window
            .document()
            .ok_or_else(|| JsValue::from_str("document is unavailable"))?;
        let canvas = document
            .get_element_by_id("emstudio-canvas")
            .ok_or_else(|| JsValue::from_str("missing canvas #emstudio-canvas"))?
            .dyn_into::<web_sys::HtmlCanvasElement>()?;

        WebRunner::new()
            .start(
                canvas,
                eframe::WebOptions::default(),
                Box::new(|_cc| Ok(Box::new(EmStudioApp::new_default()))),
            )
            .await
    }
}

#[cfg(not(target_arch = "wasm32"))]
pub fn host_build_marker() {
    // Intentionally empty: allows host cargo check for this crate.
}

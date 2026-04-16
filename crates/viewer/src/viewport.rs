use eframe::egui_wgpu;
use rmsh_renderer::Scene;

/// The egui_wgpu callback that bridges egui and the rcad-render renderer.
pub struct ViewportCallback;

impl egui_wgpu::CallbackTrait for ViewportCallback {
    fn prepare(
        &self,
        _device: &wgpu::Device,
        queue: &wgpu::Queue,
        screen_descriptor: &egui_wgpu::ScreenDescriptor,
        _encoder: &mut wgpu::CommandEncoder,
        callback_resources: &mut egui_wgpu::CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        if let Some(scene) = callback_resources.get::<Scene>() {
            scene.update_uniforms(
                queue,
                screen_descriptor.size_in_pixels[0],
                screen_descriptor.size_in_pixels[1],
            );
        }
        Vec::new()
    }

    fn paint(
        &self,
        info: egui::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        callback_resources: &egui_wgpu::CallbackResources,
    ) {
        let Some(scene) = callback_resources.get::<Scene>() else {
            return;
        };

        let viewport = info.viewport_in_pixels();
        let width = viewport.width_px as u32;
        let height = viewport.height_px as u32;

        if width == 0 || height == 0 {
            return;
        }

        render_pass.set_viewport(
            viewport.left_px as f32,
            viewport.top_px as f32,
            width as f32,
            height as f32,
            0.0,
            1.0,
        );

        // Delegate all rendering to rcad-render's WgpuRenderer
        scene.draw_in_render_pass(render_pass);

        // Draw axis gizmo in bottom-right corner
        scene.draw_axis_gizmo_in_render_pass(
            render_pass,
            [viewport.left_px as u32, viewport.top_px as u32],
            [width, height],
        );
    }
}

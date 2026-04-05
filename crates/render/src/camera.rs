use glam::{Mat4, Vec3};

pub struct OrbitCamera {
    pub target: Vec3,
    pub distance: f32,
    pub azimuth: f32,
    pub elevation: f32,
    pub fov_y: f32,
    pub near: f32,
    pub far: f32,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self {
            target: Vec3::ZERO,
            distance: 3.0,
            azimuth: 0.6,
            elevation: 0.4,
            fov_y: std::f32::consts::FRAC_PI_4,
            near: 0.01,
            far: 100.0,
        }
    }
}

impl OrbitCamera {
    pub fn rotate(&mut self, dx: f32, dy: f32) {
        self.azimuth += dx * 0.005;
        self.elevation = (self.elevation + dy * 0.005).clamp(
            -std::f32::consts::FRAC_PI_2 + 0.01,
            std::f32::consts::FRAC_PI_2 - 0.01,
        );
    }

    pub fn zoom(&mut self, delta: f32) {
        self.distance *= (1.0 - delta * 0.001).clamp(0.1, 10.0);
        self.distance = self.distance.clamp(0.1, 100.0);
    }

    pub fn eye_position(&self) -> Vec3 {
        let x = self.distance * self.elevation.cos() * self.azimuth.sin();
        let y = self.distance * self.elevation.sin();
        let z = self.distance * self.elevation.cos() * self.azimuth.cos();
        self.target + Vec3::new(x, y, z)
    }

    pub fn view_matrix(&self) -> Mat4 {
        Mat4::look_at_rh(self.eye_position(), self.target, Vec3::Y)
    }

    pub fn projection_matrix(&self, aspect: f32) -> Mat4 {
        Mat4::perspective_rh(self.fov_y, aspect, self.near, self.far)
    }

    pub fn view_projection(&self, aspect: f32) -> Mat4 {
        self.projection_matrix(aspect) * self.view_matrix()
    }

    /// Pan the camera target in screen-space.
    pub fn pan(&mut self, dx: f32, dy: f32) {
        let eye = self.eye_position();
        let forward = (self.target - eye).normalize();
        let right = forward.cross(Vec3::Y).normalize();
        let up = right.cross(forward).normalize();
        let scale = self.distance * 0.002;
        self.target += right * (-dx * scale) + up * (dy * scale);
    }

    /// Convert screen coordinates to a ray in world space (origin, direction).
    pub fn screen_to_ray(
        &self,
        x: f32,
        y: f32,
        viewport: [f32; 2],
        aspect: f32,
    ) -> ([f32; 3], [f32; 3]) {
        let inv_vp = self.view_projection(aspect).inverse();
        // Convert from viewport coords [0, width] x [0, height] to NDC [-1, 1]
        let ndc_x = (x / viewport[0]) * 2.0 - 1.0;
        let ndc_y = 1.0 - (y / viewport[1]) * 2.0; // flip Y

        let near_ndc = glam::Vec4::new(ndc_x, ndc_y, 0.0, 1.0);
        let far_ndc = glam::Vec4::new(ndc_x, ndc_y, 1.0, 1.0);

        let near_world = inv_vp * near_ndc;
        let far_world = inv_vp * far_ndc;

        let near = near_world.truncate() / near_world.w;
        let far = far_world.truncate() / far_world.w;

        let dir = (far - near).normalize();
        (near.into(), dir.into())
    }

    pub fn set_preset(&mut self, preset: ViewPreset) {
        match preset {
            ViewPreset::Front => {
                self.azimuth = 0.0;
                self.elevation = 0.0;
            }
            ViewPreset::Back => {
                self.azimuth = std::f32::consts::PI;
                self.elevation = 0.0;
            }
            ViewPreset::Left => {
                self.azimuth = -std::f32::consts::FRAC_PI_2;
                self.elevation = 0.0;
            }
            ViewPreset::Right => {
                self.azimuth = std::f32::consts::FRAC_PI_2;
                self.elevation = 0.0;
            }
            ViewPreset::Top => {
                self.azimuth = 0.0;
                self.elevation = std::f32::consts::FRAC_PI_2 - 0.01;
            }
            ViewPreset::Bottom => {
                self.azimuth = 0.0;
                self.elevation = -(std::f32::consts::FRAC_PI_2 - 0.01);
            }
            ViewPreset::Iso => {
                self.azimuth = 0.6;
                self.elevation = 0.4;
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViewPreset {
    Front,
    Back,
    Left,
    Right,
    Top,
    Bottom,
    Iso,
}

#[cfg(not(target_arch = "wasm32"))]
use eframe::NativeOptions;

#[cfg(not(target_arch = "wasm32"))]
use emstudio_app::App;
#[cfg(not(target_arch = "wasm32"))]
use emstudio_domain::Edition;
#[cfg(not(target_arch = "wasm32"))]
use emstudio_infra::RunMode;

#[cfg(not(target_arch = "wasm32"))]
fn parse_args() -> (RunMode, Edition) {
    let mut args = std::env::args().skip(1);
    let mut mode = RunMode::Standalone;
    let mut edition = Edition::Professional;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--mode" | "-m" => {
                if let Some(value) = args.next() {
                    mode = match value.as_str() {
                        "standalone" => RunMode::Standalone,
                        "cloud" => RunMode::Cloud,
                        "local-first" => RunMode::LocalFirst,
                        _ => {
                            eprintln!(
                                "[emstudio] unknown mode '{value}', fallback to standalone"
                            );
                            RunMode::Standalone
                        }
                    };
                } else {
                    eprintln!("[emstudio] missing mode value after {arg}, fallback to standalone");
                }
            }
            "--edition" | "-e" => {
                if let Some(value) = args.next() {
                    edition = match value.as_str() {
                        "basic" => Edition::Basic,
                        "professional" | "pro" => Edition::Professional,
                        "enterprise" => Edition::Enterprise,
                        _ => {
                            eprintln!(
                                "[emstudio] unknown edition '{value}', fallback to professional"
                            );
                            Edition::Professional
                        }
                    };
                } else {
                    eprintln!("[emstudio] missing edition value after {arg}, fallback to professional");
                }
            }
            _ => {}
        }
    }

    (mode, edition)
}

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result<()> {
    let (mode, edition) = parse_args();
    let options = NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };
    eframe::run_native(
        "emstudio",
        options,
        Box::new(move |cc| Ok(Box::new(App::new(mode, edition, cc)))),
    )
}

#[cfg(target_arch = "wasm32")]
fn main() {}

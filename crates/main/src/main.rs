#[cfg(not(target_arch = "wasm32"))]
use eframe::NativeOptions;

#[cfg(not(target_arch = "wasm32"))]
use emstudio_app::App;
#[cfg(not(target_arch = "wasm32"))]
use emstudio_infra::RunMode;

#[cfg(not(target_arch = "wasm32"))]
fn parse_run_mode() -> RunMode {
    let mut args = std::env::args().skip(1);
    let mut mode = RunMode::Standalone;

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
            _ => {}
        }
    }

    mode
}

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result<()> {
    let mode = parse_run_mode();
    let options = NativeOptions::default();
    eframe::run_native(
        "emstudio",
        options,
        Box::new(|_cc| Ok(Box::new(App::new(mode)))),
    )
}

#[cfg(target_arch = "wasm32")]
fn main() {}

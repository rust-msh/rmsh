//! STEP → MSH conversion tool.
//!
//! Reads a STEP file, meshes it with the specified algorithm and element size,
//! and writes the result as a Gmsh MSH file.
//!
//! Usage:
//!   cargo run -p step_to_msh -- input.step output.msh --size 0.5
//!   cargo run -p step_to_msh -- input.step output.msh --size 0.3 --algo hxt
//!   cargo run -p step_to_msh -- input.step output.msh --size 0.2 --algo frontal

use std::path::PathBuf;

use clap::Parser;
use rmsh_algo::{Delaunay3D, Frontal3D, Hxt3D, MmgRemesh, MeshParams, Mesher3D};
use rmsh_io::{load_step_from_path, save_msh_v2_to_path, save_msh_v4_to_path};

#[derive(Parser)]
#[command(name = "step_to_msh", version, about = "STEP → MSH mesh converter")]
struct Args {
    /// Input STEP file path
    input: PathBuf,

    /// Output MSH file path
    output: PathBuf,

    /// Target element edge length
    #[arg(short, long, default_value = "0.5")]
    size: f64,

    /// Meshing algorithm: delaunay, frontal, hxt, mmg
    #[arg(short, long, default_value = "hxt")]
    algo: String,

    /// MSH format version: 2 or 4
    #[arg(long, default_value = "2")]
    msh_version: u8,
}

fn main() {
    let args = Args::parse();

    // 1. Load STEP → triangulated surface
    eprintln!("Reading STEP: {}", args.input.display());
    let surface: rmsh_model::Mesh = match load_step_from_path(&args.input) {
        Ok(m) => m,
        Err(e) => { eprintln!("ERROR: failed to read STEP: {e}"); std::process::exit(1); }
    };
    eprintln!("  surface: {} nodes, {} triangles", surface.node_count(), surface.element_count());

    // 2. Select 3D mesher
    let params = MeshParams::with_size(args.size);
    let mesher: Box<dyn Mesher3D> = match args.algo.as_str() {
        "delaunay" => Box::new(Delaunay3D::default()),
        "frontal"  => Box::new(Frontal3D::default()),
        "hxt"      => Box::new(Hxt3D::default()),
        "mmg"      => Box::new(MmgRemesh::default()),
        other => { eprintln!("ERROR: unknown algorithm '{other}' (delaunay/frontal/hxt/mmg)");
                   std::process::exit(1); },
    };

    // 3. Generate tetrahedral mesh
    eprintln!("Meshing (algo={}, size={})...", args.algo, args.size);
    let mesh = match mesher.mesh_3d(&surface, &params) {
        Ok(m) => m,
        Err(e) => { eprintln!("ERROR: meshing failed: {e}"); std::process::exit(1); }
    };
    eprintln!("  result: {} nodes, {} tets", mesh.node_count(),
              mesh.elements_by_dimension(3).len());

    // 4. Write MSH output
    eprintln!("Writing MSH v{}: {}", args.msh_version, args.output.display());
    let result = match args.msh_version {
        2 => save_msh_v2_to_path(&args.output, &mesh),
        4 => save_msh_v4_to_path(&args.output, &mesh),
        v => { eprintln!("ERROR: unsupported MSH version {v}"); std::process::exit(1); }
    };
    match result {
        Ok(()) => eprintln!("Done."),
        Err(e) => { eprintln!("ERROR: failed to write MSH: {e}"); std::process::exit(1); }
    }
}

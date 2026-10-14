mod msh;
mod msh_partitioned;
pub mod vtk;
pub use vtk::{read_vtk_legacy, read_vtu, write_vtk_legacy, write_vtu};
mod msh_periodic;

pub use msh::{
    MshError, load_msh_from_bytes, load_msh_from_path, parse_msh, save_msh_v2_binary_to_path,
    save_msh_v2_to_path, save_msh_v4_binary_to_path, save_msh_v4_to_path, write_msh_v2,
    write_msh_v2_binary, write_msh_v4, write_msh_v4_binary,
};
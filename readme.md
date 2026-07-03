# rmsh

A Rust mesh generation and optimization framework with CAD kernel integration.

## Features

- 2D/3D mesh generation (Delaunay, frontal, paving)
- Mesh optimization (Laplacian smoothing, quality optimization)
- CAD boolean operations (union, intersection, difference) via rcad2
- STEP and Gmsh MSH file I/O
- Interactive 3D viewer (desktop via egui + wgpu)
- Python bindings (PyO3)

## Project Structure

```
rmsh/
├── crates/
│   ├── algo/        # Mesh generation & optimization algorithms
│   ├── model/       # Mesh data structures (Node, Element, GModel)
│   ├── geo/         # Geometry processing & topology classification
│   ├── io/          # File I/O (STEP, Gmsh MSH v2/v4)
│   ├── renderer/    # WebGPU rendering pipeline
│   ├── viewer/      # Desktop viewer application
│   └── py/          # Python bindings
├── vendor/
│   └── rcad2/       # CAD kernel submodule (BRep, boolean ops, STEP)
│       ├── libs/rcad-kernel/       # BRep data structures
│       ├── libs/rcad-modeling/     # Primitive shape builders
│       ├── libs/rcad-algorithms/   # Boolean operations
│       ├── libs/rcad-step/         # STEP reader/writer
│       └── libs/rcad-render/       # GPU tessellation & rendering
└── testdata/
```

## Run viewer as desktop app

```
cargo run -p rmsh-viewer
```

## Run viewer as web app

```
trunk serve
```

## Run tests

```sh
cargo test --workspace
```

## STEP → MSH command-line tool

```sh
# Convert a STEP file to a tetrahedral mesh (MSH v2)
cargo run -p step_to_msh -- input.step output.msh --size 0.5

# Specify algorithm and MSH version
cargo run -p step_to_msh -- input.step output.msh --size 0.3 --algo hxt --msh-version 4

# Available algorithms: delaunay, frontal, hxt (default), mmg
```

Output can be opened in Gmsh or rmsh-viewer for visual inspection.

## Python bindings

```sh
cd crates/py
maturin develop
python -c "import rmsh; rmsh.initialize()"
```
from pathlib import Path

import rmsh


def export_solids_only_one_step(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "solids_only.step"

    rmsh.initialize()
    try:
        rmsh.clear()

        # Box
        rmsh.model.occ.addBox(0.0, 0.0, 0.0, 2.0, 1.2, 1.0)

        # Sphere
        rmsh.model.occ.addSphere(3.8, 0.8, 0.8, 0.8)

        # Cone
        rmsh.model.occ.addCone(6.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.9, 0.2)

        # Torus
        rmsh.model.occ.addTorus(8.8, 0.8, 0.8, 0.0, 0.0, 1.0, 0.9, 0.25)

        # One-shot export
        rmsh.write(str(out_file))
    finally:
        rmsh.finalize()

    print(f"Exported: {out_file}")


if __name__ == "__main__":
    output_dir = Path(__file__).resolve().parent / "out_solids_only_step"
    export_solids_only_one_step(output_dir)

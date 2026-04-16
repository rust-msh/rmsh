"""API upgrade smoke test for rmsh Python bindings.

Covers newly implemented logger/model/physical/plugin/mesh utility APIs.
"""

import rmsh


def main() -> None:
    rmsh.initialize()
    rmsh.logger.start()

    rmsh.model.add("api_upgrade_demo")
    assert rmsh.model.getCurrent() == "api_upgrade_demo"

    # Build a simple CAD surface and attach model metadata.
    rmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 2.0, 1.0)
    rmsh.model.setEntityName(2, 1, "demo_surface")
    assert rmsh.model.getEntityName(2, 1) == "demo_surface"

    pg = rmsh.model.addPhysicalGroup(2, [1], -1, "demo_group")
    assert (2, pg) in rmsh.model.getPhysicalGroups(-1)
    assert rmsh.model.getPhysicalName(2, pg) == "demo_group"
    rmsh.model.setPhysicalName(2, pg, "demo_group_renamed")
    assert rmsh.model.getPhysicalName(2, pg) == "demo_group_renamed"

    # Geometry-first export path (CAD/BRep -> STEP).
    rmsh.write("api_upgrade_geom.step")

    rmsh.model.mesh.generate(2)

    # Exercise mesh utilities.
    refined = rmsh.model.mesh.refine()
    recombined = rmsh.model.mesh.recombine(2, -1, 45.0)
    print(f"refined triangles: {refined}, recombined quads: {recombined}")

    rmsh.model.mesh.setOrder(2)

    # Plugin API dispatch.
    rmsh.plugin.setNumber("Smooth", "niter", 5)
    rmsh.plugin.setString("Smooth", "method", "Laplace")
    rmsh.plugin.run("Smooth")

    rmsh.write("api_upgrade_mesh.msh")

    logs = rmsh.logger.get()
    print(f"logger entries: {len(logs)}")
    if logs:
        print("last log:", logs[-1])

    rmsh.logger.stop()
    rmsh.finalize()
    print("api upgrade example passed")


if __name__ == "__main__":
    main()

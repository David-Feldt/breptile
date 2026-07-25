import os

import pytest


@pytest.fixture(scope="session")
def fixtures(tmp_path_factory):
    """Generate STL fixtures from known build123d models."""
    from build123d import Box, Cylinder, Mesher, Part, Pos, Sphere, Torus

    d = tmp_path_factory.mktemp("fixtures")

    def write(name, part, deflection=0.05):
        path = os.path.join(d, name)
        m = Mesher()
        m.add_shape(part, linear_deflection=deflection, angular_deflection=0.3)
        m.write(path)
        return path

    out = {}
    out["cube"] = write("cube.stl", Box(20, 20, 20))
    out["plate_hole"] = write(
        "plate_hole.stl", Box(40, 30, 8) - Cylinder(5, 8)
    )
    out["torus"] = write("torus.stl", Torus(20, 5))
    out["sphere"] = write("sphere.stl", Sphere(10), deflection=0.02)
    out["ball_pocket"] = write(
        "ball_pocket.stl", Box(30, 30, 10) - Pos(0, 0, 5) * Sphere(6), deflection=0.02
    )
    return out

import numpy as np
import pytest

from stl2step.mesh import MeshError, load_mesh
from stl2step.pipeline import convert


def _cyl_face_count(step_path):
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_SurfaceType
    from OCP.TopAbs import TopAbs_ShapeEnum
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from build123d import import_step

    shape = import_step(step_path).wrapped
    n = 0
    exp = TopExp_Explorer(shape, TopAbs_ShapeEnum.TopAbs_FACE)
    while exp.More():
        surf = BRepAdaptor_Surface(TopoDS.Face_s(exp.Current()))
        if surf.GetType() == GeomAbs_SurfaceType.GeomAbs_Cylinder:
            n += 1
        exp.Next()
    return n


def test_cube_tessellated(fixtures, tmp_path):
    out = str(tmp_path / "cube.step")
    report = convert(fixtures["cube"], out, mode="tessellated")
    assert report["valid_brep"]
    # unify should collapse the 12 triangles to 6 planar faces
    assert report["faces"] == 6


def test_cube_fit(fixtures, tmp_path):
    out = str(tmp_path / "cube_fit.step")
    report = convert(fixtures["cube"], out, mode="fit")
    assert report["valid_brep"]
    assert report["faces"] == 6


def test_plate_hole_fit(fixtures, tmp_path):
    out = str(tmp_path / "plate.step")
    report = convert(fixtures["plate_hole"], out, mode="fit")
    assert report["valid_brep"]
    # 6 planar box faces + 1 cylindrical hole face
    assert report["faces"] == 7
    assert _cyl_face_count(out) == 1


def test_plate_hole_deviation(fixtures, tmp_path):
    from stl2step.verify import deviation

    out = str(tmp_path / "plate_dev.step")
    convert(fixtures["plate_hole"], out, mode="fit")
    mesh = load_mesh(fixtures["plate_hole"])
    dev = deviation(out, mesh, samples=5000)
    # STEP surface may legitimately deviate from facets by the chord error;
    # bound by tessellation deflection used for the fixture.
    assert dev["max"] < 0.2


def test_torus_freeform_fallback(fixtures, tmp_path):
    out = str(tmp_path / "torus.step")
    report = convert(fixtures["torus"], out, mode="fit")
    assert report["faces"] > 0  # tessellated fallback, but a valid STEP emerges


def test_full_sphere(fixtures, tmp_path):
    out = str(tmp_path / "sphere.step")
    report = convert(fixtures["sphere"], out, mode="fit")
    assert report["valid_brep"]
    assert report["faces"] == 1


def test_spherical_pocket(fixtures, tmp_path):
    out = str(tmp_path / "pocket.step")
    report = convert(fixtures["ball_pocket"], out, mode="fit")
    assert report["valid_brep"]
    assert report["regions"]["sphere"] == 1
    assert report["regions"]["fallback"] == 0


def test_non_watertight_errors(tmp_path):
    import trimesh

    # A single open triangle sheet.
    m = trimesh.Trimesh(vertices=[[0, 0, 0], [1, 0, 0], [0, 1, 0]], faces=[[0, 1, 2]])
    p = str(tmp_path / "open.stl")
    m.export(p)
    with pytest.raises(MeshError):
        load_mesh(p)
    load_mesh(p, force=True)  # allowed with force

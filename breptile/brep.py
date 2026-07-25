"""Build OCC BREP shapes from mesh data via OCP."""

from __future__ import annotations

import numpy as np

from OCP.BRep import BRep_Builder
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakePolygon,
    BRepBuilderAPI_MakeSolid,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_Sewing,
)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Pnt, gp_Vec
from OCP.ShapeFix import ShapeFix_Shape, ShapeFix_Solid
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Compound, TopoDS_Shape, TopoDS_Shell


def _pnt(v) -> gp_Pnt:
    return gp_Pnt(float(v[0]), float(v[1]), float(v[2]))


def triangle_face(a, b, c) -> TopoDS_Shape | None:
    """Planar face from one triangle; None if degenerate."""
    poly = BRepBuilderAPI_MakePolygon(_pnt(a), _pnt(b), _pnt(c), True)
    if not poly.IsDone():
        return None
    face = BRepBuilderAPI_MakeFace(poly.Wire(), True)
    return face.Face() if face.IsDone() else None


def polygon_wire(points: np.ndarray) -> TopoDS_Shape | None:
    """Closed polygonal wire through ordered 3D points."""
    poly = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly.Add(_pnt(p))
    poly.Close()
    return poly.Wire() if poly.IsDone() else None


def sew(faces: list[TopoDS_Shape], tolerance: float) -> TopoDS_Shape:
    sewing = BRepBuilderAPI_Sewing(tolerance)
    for f in faces:
        sewing.Add(f)
    sewing.Perform()
    return sewing.SewedShape()


def shells_of(shape: TopoDS_Shape) -> list[TopoDS_Shell]:
    out = []
    exp = TopExp_Explorer(shape, TopAbs_ShapeEnum.TopAbs_SHELL)
    while exp.More():
        out.append(TopoDS.Shell_s(exp.Current()))
        exp.Next()
    return out


def to_solid(sewed: TopoDS_Shape) -> TopoDS_Shape:
    """Promote sewn shell(s) to a solid; falls back to the shell on failure."""
    shells = shells_of(sewed)
    if not shells:
        return sewed
    maker = BRepBuilderAPI_MakeSolid()
    for sh in shells:
        maker.Add(sh)
    if not maker.IsDone():
        return sewed
    solid = maker.Solid()
    fix = ShapeFix_Solid(solid)
    fix.Perform()
    fixed = fix.Solid()
    # Guard against inside-out solids: volume must be positive.
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(fixed, props)
    if props.Mass() < 0:
        fixed.Reverse()
    return fixed


def unify(shape: TopoDS_Shape, linear_tol: float, angular_tol: float) -> TopoDS_Shape:
    u = ShapeUpgrade_UnifySameDomain(shape, True, True, False)
    u.SetLinearTolerance(linear_tol)
    u.SetAngularTolerance(angular_tol)
    u.Build()
    return u.Shape()


def heal(shape: TopoDS_Shape) -> TopoDS_Shape:
    fix = ShapeFix_Shape(shape)
    fix.Perform()
    return fix.Shape()


def is_valid(shape: TopoDS_Shape) -> bool:
    return BRepCheck_Analyzer(shape).IsValid()


def count_faces(shape: TopoDS_Shape) -> int:
    n = 0
    exp = TopExp_Explorer(shape, TopAbs_ShapeEnum.TopAbs_FACE)
    while exp.More():
        n += 1
        exp.Next()
    return n


def compound_of(shapes: list[TopoDS_Shape]) -> TopoDS_Compound:
    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    for s in shapes:
        builder.Add(comp, s)
    return comp


def surface_patch_face(
    surface, loops_pts: list[np.ndarray], tolerance: float
) -> TopoDS_Shape | None:
    """Face on an analytic ``Geom_Surface`` trimmed by closed polyline loops.

    Loops are given as ordered 3D point arrays lying on (or near) the surface;
    pcurves are obtained by projection, so the caller must orient the surface's
    parametric seam away from the patch.
    """
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.ShapeFix import ShapeFix_Edge, ShapeFix_Face
    from OCP.TopAbs import TopAbs_ShapeEnum
    from OCP.TopExp import TopExp_Explorer

    base = BRepBuilderAPI_MakeFace(surface, tolerance)
    if not base.IsDone():
        return None
    base_face = base.Face()

    wires = []
    for pts in loops_pts:
        wire = polygon_wire(pts)
        if wire is None:
            return None
        fixer = ShapeFix_Edge()
        exp = TopExp_Explorer(wire, TopAbs_ShapeEnum.TopAbs_EDGE)
        while exp.More():
            fixer.FixAddPCurve(TopoDS.Edge_s(exp.Current()), base_face, False, tolerance)
            exp.Next()
        wires.append(wire)

    maker = BRepBuilderAPI_MakeFace(surface, wires[0], False)
    if not maker.IsDone():
        return None
    for w in wires[1:]:
        maker.Add(w)
    if not maker.IsDone():
        return None
    fix = ShapeFix_Face(maker.Face())
    fix.Perform()
    face = fix.Face()
    return face if BRepCheck_Analyzer(face).IsValid() else None


def mesh_to_tessellated_solid(mesh, tolerance: float) -> TopoDS_Shape:
    """One planar face per triangle → sew → solid. The guaranteed-success path."""
    verts = np.asarray(mesh.vertices)
    faces = []
    for tri in mesh.faces:
        f = triangle_face(verts[tri[0]], verts[tri[1]], verts[tri[2]])
        if f is not None:
            faces.append(f)
    if not faces:
        raise RuntimeError("no valid faces could be built from mesh")
    return to_solid(sew(faces, tolerance))

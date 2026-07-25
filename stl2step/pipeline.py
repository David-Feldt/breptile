"""Conversion pipeline: mesh → regions → BREP → STEP."""

from __future__ import annotations

import numpy as np

from stl2step import brep
from stl2step.mesh import load_mesh
from stl2step.segment import Region, boundary_loops, segment


def _default_tol(mesh) -> float:
    return float(np.linalg.norm(mesh.bounding_box.extents) * 1e-4)


def _planar_region_face(mesh, region: Region, circle_subs: dict):
    """Single planar face for a region, holes included. Loops whose vertices all
    belong to a replaced cylinder rim become true circles (from circle_subs)."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.ShapeFix import ShapeFix_Face
    from OCP.gp import gp_Ax3, gp_Dir, gp_Pln, gp_Pnt

    loops = boundary_loops(mesh, region.faces)
    if not loops:
        return None
    o, n = region.params["origin"], region.params["normal"]
    pln = gp_Pln(
        gp_Ax3(gp_Pnt(*map(float, o)), gp_Dir(*map(float, n)))
    )

    wires = []
    areas = []
    for loop in loops:
        wire = None
        key = frozenset(loop.tolist())
        if key in circle_subs:
            wire = circle_subs[key]
        if wire is None:
            wire = brep.polygon_wire(mesh.vertices[loop])
        if wire is None:
            return None
        pts = mesh.vertices[loop]
        # Loop area via shoelace on the plane (for outer-loop detection).
        rel = pts - o
        e1 = np.cross(n, [1.0, 0, 0])
        if np.linalg.norm(e1) < 1e-6:
            e1 = np.cross(n, [0, 1.0, 0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(n, e1)
        uv = np.column_stack([rel @ e1, rel @ e2])
        area = abs(
            0.5 * np.sum(uv[:, 0] * np.roll(uv[:, 1], -1) - np.roll(uv[:, 0], -1) * uv[:, 1])
        )
        wires.append(wire)
        areas.append(area)

    order = np.argsort(areas)[::-1]
    maker = BRepBuilderAPI_MakeFace(pln, wires[order[0]], True)
    if not maker.IsDone():
        return None
    for i in order[1:]:
        maker.Add(wires[i])
    if not maker.IsDone():
        return None
    face = maker.Face()
    fix = ShapeFix_Face(face)
    fix.FixOrientation()
    return fix.Face()


def _cylinder_region_faces(mesh, region: Region, tol: float):
    """Analytic cylindrical face for a full-wrap cylinder region.

    Returns (faces, circle_subs) or (None, {}) when the region doesn't wrap the
    full circumference (partial fillets etc. stay tessellated in v1).
    circle_subs maps rim-loop vertex sets → circular wire, so adjacent planar
    faces can share the exact same rim geometry and sewing closes the gap.
    """
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
    from OCP.gp import gp_Ax2, gp_Circ, gp_Dir, gp_Pnt, gp_Vec

    axis = region.params["axis"]
    center = region.params["center"]
    r = region.params["radius"]

    loops = boundary_loops(mesh, region.faces)

    # Angular coverage of region vertices around the axis (full wrap check).
    pts = mesh.vertices[np.unique(mesh.faces[region.faces])]
    rel = pts - center
    e1 = np.cross(axis, [1.0, 0, 0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.cross(axis, [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    ang = np.arctan2(rel @ e2, rel @ e1)
    hist = np.histogram(ang, bins=24, range=(-np.pi, np.pi))[0]
    if (hist == 0).any() or len(loops) != 2:
        return _partial_cylinder_faces(mesh, region, loops, ang, e1, e2, tol)

    circle_subs = {}
    rim_wires = []
    for loop in loops:
        lp = mesh.vertices[loop]
        # Rim must be circular: constant height and radius within tol.
        t = (lp - center) @ axis
        rad = np.linalg.norm((lp - center) - np.outer(t, axis), axis=1)
        if np.ptp(t) > tol * 4 or np.abs(rad - r).max() > tol * 4:
            return None, {}
        h = t.mean()
        c = center + h * axis
        circ = gp_Circ(gp_Ax2(gp_Pnt(*map(float, c)), gp_Dir(*map(float, axis))), float(r))
        edge = BRepBuilderAPI_MakeEdge(circ).Edge()
        wire = BRepBuilderAPI_MakeWire(edge).Wire()
        circle_subs[frozenset(loop.tolist())] = wire
        rim_wires.append((h, wire))

    from OCP.gp import gp_Ax3, gp_Cylinder

    h0 = min(rw[0] for rw in rim_wires)
    h1 = max(rw[0] for rw in rim_wires)
    if h1 - h0 < tol:
        return None, {}
    base = center + h0 * axis
    cyl = gp_Cylinder(gp_Ax3(gp_Pnt(*map(float, base)), gp_Dir(*map(float, axis))), float(r))
    maker = BRepBuilderAPI_MakeFace(cyl, 0.0, 2 * np.pi, 0.0, float(h1 - h0))
    if not maker.IsDone():
        return None, {}
    return [maker.Face()], circle_subs


def _partial_cylinder_faces(mesh, region, loops, ang, e1, e2, tol):
    """Trimmed patch on a cylinder that doesn't wrap the full circumference.
    The parametric seam (u=0) is rotated into the angular gap so pcurve
    projection never crosses it."""
    from OCP.Geom import Geom_CylindricalSurface
    from OCP.gp import gp_Ax3, gp_Dir, gp_Pnt

    if not loops:
        return None, {}
    axis, center, r = region.params["axis"], region.params["center"], region.params["radius"]

    # Find the largest angular gap; point X of the frame at its middle.
    hist, edges = np.histogram(ang, bins=48, range=(-np.pi, np.pi))
    empty = np.nonzero(hist == 0)[0]
    if len(empty) == 0:
        return None, {}  # full wrap but wrong loop count — leave tessellated
    seam_ang = (edges[empty[0]] + edges[empty[0] + 1]) / 2
    xdir = np.cos(seam_ang) * e1 + np.sin(seam_ang) * e2

    frame = gp_Ax3(gp_Pnt(*map(float, center)), gp_Dir(*map(float, axis)),
                   gp_Dir(*map(float, xdir)))
    surf = Geom_CylindricalSurface(frame, float(r))
    face = brep.surface_patch_face(surf, [mesh.vertices[lp] for lp in loops], tol * 4)
    return ([face], {}) if face is not None else (None, {})


def _circular_loop(pts, tol):
    """If a loop is a planar circle, return (axis, center, radius); else None."""
    centroid = pts.mean(axis=0)
    _, s, vt = np.linalg.svd(pts - centroid, full_matrices=False)
    normal = vt[-1]
    if np.abs((pts - centroid) @ normal).max() > tol * 4:
        return None
    radii = np.linalg.norm((pts - centroid) - np.outer((pts - centroid) @ normal, normal), axis=1)
    if np.ptp(radii) > tol * 8:
        return None
    return normal, centroid, float(radii.mean())


def _sphere_region_faces(mesh, region, tol):
    """Full spheres → one face; bands/countersink zones wrapping an axis →
    v-trimmed zone with true circular rims (returned via circle_subs, like
    full cylinders); non-wrapping caps → projected patch.

    Returns (faces | None, circle_subs).
    """
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
    from OCP.Geom import Geom_SphericalSurface
    from OCP.gp import gp_Ax2, gp_Ax3, gp_Circ, gp_Dir, gp_Pnt, gp_Sphere

    center, r = region.params["center"], region.params["radius"]
    loops = boundary_loops(mesh, region.faces)

    if not loops:
        frame = gp_Ax3(gp_Pnt(*map(float, center)), gp_Dir(0, 0, 1))
        maker = BRepBuilderAPI_MakeFace(gp_Sphere(frame, float(r)))
        return ([maker.Face()], {}) if maker.IsDone() else (None, {})

    # Band/zone: every loop is a circle about a common axis.
    circles = [_circular_loop(mesh.vertices[lp], tol) for lp in loops]
    if all(c is not None for c in circles) and len(circles) in (1, 2):
        axis = circles[0][0]
        if all(abs(np.dot(axis, c[0])) > 0.999 for c in circles):
            # Orient axis away from center through the patch.
            pts = mesh.vertices[np.unique(mesh.faces[region.faces])]
            d = pts.mean(axis=0) - center
            if np.dot(axis, d) < 0:
                axis = -axis
            heights = sorted(float(np.clip((c[1] - center) @ axis, -r, r)) for c in circles)
            v = [float(np.arcsin(np.clip(h / r, -1, 1))) for h in heights]
            # Single rim → cap from the rim up to the pole the patch points at
            # (the frame axis was oriented through the patch centroid above).
            v0, v1 = (v[0], v[1]) if len(v) == 2 else (v[0], np.pi / 2)
            if v1 - v0 < 1e-6:
                return None, {}
            # The mesh must actually cover the u×v extent we are about to
            # claim, or the face would invent surface (e.g. a rim loop on a
            # partial patch would get promoted to a full cap).
            rel = pts - center
            e1 = np.cross(axis, [1.0, 0, 0])
            if np.linalg.norm(e1) < 1e-6:
                e1 = np.cross(axis, [0, 1.0, 0])
            e1 /= np.linalg.norm(e1)
            e2 = np.cross(axis, e1)
            u_pts = np.arctan2(rel @ e2, rel @ e1)
            v_pts = np.arcsin(np.clip((rel @ axis) / r, -1, 1))
            u_hist = np.histogram(u_pts, bins=16, range=(-np.pi, np.pi))[0]
            margin = max((v1 - v0) * 0.35, 0.15)
            covered = (u_hist > 0).all() and v_pts.min() < v0 + margin and v_pts.max() > v1 - margin
            if not covered:
                return None, {}
            frame = gp_Ax3(gp_Pnt(*map(float, center)), gp_Dir(*map(float, axis)))
            maker = BRepBuilderAPI_MakeFace(gp_Sphere(frame, float(r)), 0.0, 2 * np.pi, v0, v1)
            if maker.IsDone():
                subs = {}
                for lp, (n_, c_, rad_) in zip(loops, circles):
                    h = float(np.clip((c_ - center) @ axis, -r, r))
                    rim_c = center + h * axis
                    rim_r = float(np.sqrt(max(r * r - h * h, 1e-12)))
                    circ = gp_Circ(gp_Ax2(gp_Pnt(*map(float, rim_c)), gp_Dir(*map(float, axis))), rim_r)
                    wire = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire()
                    subs[frozenset(lp.tolist())] = wire
                return [maker.Face()], subs

    # Non-wrapping patch: put its centroid on the equator, seam and poles away.
    d = mesh.vertices[np.unique(mesh.faces[region.faces])].mean(axis=0) - center
    n = np.linalg.norm(d)
    if n < 1e-9:
        return None, {}
    d /= n
    helper = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
    zdir = np.cross(d, helper)
    zdir /= np.linalg.norm(zdir)
    frame = gp_Ax3(gp_Pnt(*map(float, center)), gp_Dir(*map(float, zdir)),
                   gp_Dir(*map(float, -d)))
    surf = Geom_SphericalSurface(frame, float(r))
    face = brep.surface_patch_face(surf, [mesh.vertices[lp] for lp in loops], tol * 4)
    return ([face], {}) if face is not None else (None, {})


def _tessellated_faces(mesh, face_indices):
    verts = mesh.vertices
    out = []
    for tri in mesh.faces[face_indices]:
        f = brep.triangle_face(verts[tri[0]], verts[tri[1]], verts[tri[2]])
        if f is not None:
            out.append(f)
    return out


def rebuild(mesh, regions: list[Region], tol: float, report: dict) -> "brep.TopoDS_Shape":
    """Fit-mode rebuild: analytic faces where fits are confident, facets elsewhere."""
    faces = []
    circle_subs: dict = {}
    stats = {"plane": 0, "cylinder": 0, "sphere": 0, "freeform_triangles": 0, "fallback": 0}

    cylinder_regions = [r for r in regions if r.kind == "cylinder"]
    sphere_regions = [r for r in regions if r.kind == "sphere"]
    planar_regions = [r for r in regions if r.kind == "plane"]
    other = [r for r in regions if r.kind not in ("cylinder", "sphere", "plane")]

    for r in cylinder_regions:
        built, subs = _cylinder_region_faces(mesh, r, tol)
        if built:
            faces += built
            circle_subs.update(subs)
            stats["cylinder"] += 1
        else:
            faces += _tessellated_faces(mesh, r.faces)
            stats["fallback"] += 1
            r.kind = "freeform"

    for r in sphere_regions:
        built, subs = _sphere_region_faces(mesh, r, tol)
        if built:
            faces += built
            circle_subs.update(subs)
            stats["sphere"] += 1
        else:
            faces += _tessellated_faces(mesh, r.faces)
            stats["freeform_triangles"] += len(r.faces)
            stats["fallback"] += 1

    for r in planar_regions:
        f = _planar_region_face(mesh, r, circle_subs)
        if f is not None:
            faces.append(f)
            stats["plane"] += 1
        else:
            faces += _tessellated_faces(mesh, r.faces)
            stats["fallback"] += 1

    for r in other:
        faces += _tessellated_faces(mesh, r.faces)
        stats["freeform_triangles"] += len(r.faces)

    report["regions"] = stats
    shape = brep.to_solid(brep.sew(faces, tol))
    return shape


def convert(
    input_path: str,
    output_path: str,
    mode: str = "auto",
    tol: float | None = None,
    schema: str = "AP214",
    force: bool = False,
    max_triangles: int | None = None,
) -> dict:
    """Run the pipeline. Returns a report dict."""
    from stl2step.export import write_step

    mesh = load_mesh(input_path, force=force, max_triangles=max_triangles)
    tol = tol if tol is not None else _default_tol(mesh)
    report: dict = {
        "input": input_path,
        "output": output_path,
        "mode": mode,
        "tolerance": tol,
        "triangles": int(len(mesh.faces)),
        "watertight": bool(mesh.is_watertight),
    }

    if mode == "tessellated":
        shape = brep.mesh_to_tessellated_solid(mesh, tol)
        shape = brep.unify(shape, tol, np.radians(0.5))
    else:
        regions = segment(mesh, tol)
        if mode == "prismatic":
            for r in regions:
                if r.kind in ("cylinder", "sphere"):
                    r.kind = "freeform"
        report["segmentation"] = {
            k: sum(1 for r in regions if r.kind == k)
            for k in ("plane", "cylinder", "sphere", "freeform")
        }
        report["fit_residuals"] = {
            "max_plane": max((r.residual for r in regions if r.kind == "plane"), default=0.0),
            "max_cylinder": max((r.residual for r in regions if r.kind == "cylinder"), default=0.0),
        }
        shape = rebuild(mesh, regions, tol, report)
        shape = brep.unify(shape, tol, np.radians(0.5))

    shape = brep.heal(shape)
    report["valid_brep"] = brep.is_valid(shape)
    report["faces"] = brep.count_faces(shape)
    write_step(shape, output_path, schema=schema)
    return report

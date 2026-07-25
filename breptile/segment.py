"""Mesh segmentation into planar / cylindrical / spherical / freeform regions,
with least-squares primitive fitting validated against a deviation tolerance."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import trimesh


@dataclass
class Region:
    kind: str  # "plane" | "cylinder" | "sphere" | "freeform"
    faces: np.ndarray  # triangle indices into mesh.faces
    params: dict = field(default_factory=dict)
    residual: float = 0.0  # max deviation of region vertices from fitted surface


# ---------------------------------------------------------------- fitting

def fit_plane(points: np.ndarray) -> tuple[dict, float]:
    centroid = points.mean(axis=0)
    _, _, vt = np.linalg.svd(points - centroid, full_matrices=False)
    normal = vt[-1]
    d = (points - centroid) @ normal
    return {"origin": centroid, "normal": normal}, float(np.abs(d).max())


def fit_sphere(points: np.ndarray) -> tuple[dict, float]:
    # Linear least squares: |p|^2 = 2 c·p + (r^2 - |c|^2)
    A = np.hstack([2 * points, np.ones((len(points), 1))])
    b = (points**2).sum(axis=1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    center = sol[:3]
    r = float(np.sqrt(sol[3] + center @ center))
    dev = np.abs(np.linalg.norm(points - center, axis=1) - r)
    return {"center": center, "radius": r}, float(dev.max())


def fit_cylinder(points: np.ndarray, normals: np.ndarray) -> tuple[dict, float]:
    """Axis from face normals (all ⟂ axis), circle fit in the projected plane,
    then nonlinear refinement."""
    from scipy.optimize import least_squares

    # Smallest principal direction of the normal cloud is the axis.
    _, _, vt = np.linalg.svd(normals - normals.mean(axis=0), full_matrices=False)
    axis = vt[-1]
    axis /= np.linalg.norm(axis)

    def project_basis(ax):
        e1 = np.cross(ax, [1.0, 0.0, 0.0])
        if np.linalg.norm(e1) < 1e-6:
            e1 = np.cross(ax, [0.0, 1.0, 0.0])
        e1 /= np.linalg.norm(e1)
        return e1, np.cross(ax, e1)

    centroid = points.mean(axis=0)

    def residuals(x):
        ax = x[:3] / np.linalg.norm(x[:3])
        c = x[3:6]
        r = x[6]
        rel = points - c
        radial = rel - np.outer(rel @ ax, ax)
        return np.linalg.norm(radial, axis=1) - r

    # Initial circle via Kasa fit in projection.
    e1, e2 = project_basis(axis)
    rel = points - centroid
    uv = np.column_stack([rel @ e1, rel @ e2])
    A = np.hstack([2 * uv, np.ones((len(uv), 1))])
    b = (uv**2).sum(axis=1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    c2, r0 = sol[:2], float(np.sqrt(max(sol[2] + sol[:2] @ sol[:2], 1e-12)))
    c0 = centroid + c2[0] * e1 + c2[1] * e2

    x0 = np.concatenate([axis, c0, [r0]])
    try:
        res = least_squares(residuals, x0, method="lm", max_nfev=200)
        x = res.x
    except Exception:
        x = x0
    ax = x[:3] / np.linalg.norm(x[:3])
    c, r = x[3:6], float(abs(x[6]))
    # Anchor center onto the axis near the data, record extent along axis.
    t = (points - c) @ ax
    c = c + t.mean() * ax
    dev = np.abs(residuals(np.concatenate([ax, c, [r]])))
    return {
        "axis": ax,
        "center": c,
        "radius": r,
        "half_length": float((t.max() - t.min()) / 2 + 1e-9),
    }, float(dev.max())


# ---------------------------------------------------------------- segmentation

def _smooth_components(mesh: trimesh.Trimesh, angle: float) -> list[np.ndarray]:
    """Connected components of triangles whose shared-edge dihedral deviation
    is below ``angle`` radians — one component per smooth surface patch."""
    adj = mesh.face_adjacency
    ok = mesh.face_adjacency_angles < angle
    return [
        np.asarray(c)
        for c in trimesh.graph.connected_components(
            adj[ok], min_len=1, nodes=np.arange(len(mesh.faces))
        )
    ]


def segment(mesh: trimesh.Trimesh, tol: float, smooth_angle: float = np.radians(30)) -> list[Region]:
    """Grow smooth patches, then classify each: plane → cylinder → sphere → freeform."""
    regions = []
    for comp in _smooth_components(mesh, smooth_angle):
        pts = mesh.vertices[np.unique(mesh.faces[comp])]
        normals = mesh.face_normals[comp]
        best: Region | None = None

        params, residual = fit_plane(pts)
        if residual <= tol:
            best = Region("plane", comp, params, residual)
        elif len(comp) >= 6 and len(pts) >= 8:
            cyl_params, cyl_res = fit_cylinder(pts, normals)
            if cyl_res <= tol:
                best = Region("cylinder", comp, cyl_params, cyl_res)
            else:
                sph_params, sph_res = fit_sphere(pts)
                if sph_res <= tol:
                    best = Region("sphere", comp, sph_params, sph_res)

        if best is None:
            best = Region("freeform", comp, {}, float("nan"))
        regions.append(best)
    return regions


# ---------------------------------------------------------------- boundaries

def boundary_loops(mesh: trimesh.Trimesh, face_indices: np.ndarray) -> list[np.ndarray]:
    """Ordered vertex-index loops bounding a region, following triangle winding."""
    tris = mesh.faces[face_indices]
    directed = {}
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            directed[(u, v)] = True
    nxt = {}
    for (u, v) in directed:
        if (v, u) not in directed:  # boundary half-edge
            nxt[u] = v
    loops = []
    visited = set()
    for start in list(nxt):
        if start in visited:
            continue
        loop = [start]
        visited.add(start)
        cur = nxt[start]
        while cur != start and cur in nxt and cur not in visited:
            loop.append(cur)
            visited.add(cur)
            cur = nxt[cur]
        if cur == start and len(loop) >= 3:
            loops.append(np.asarray(loop))
    return loops

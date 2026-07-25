"""Verification: deviation between an exported STEP and the source mesh."""

from __future__ import annotations

import numpy as np
import trimesh


def step_to_mesh(step_path: str, linear_deflection: float = 0.05) -> trimesh.Trimesh:
    from build123d import import_step

    shape = import_step(step_path)
    verts, faces = shape.tessellate(linear_deflection)
    return trimesh.Trimesh(
        vertices=np.asarray([(v.X, v.Y, v.Z) for v in verts]),
        faces=np.asarray(faces),
        process=False,
    )


def deviation(step_path: str, mesh: trimesh.Trimesh, samples: int = 20000) -> dict:
    """Two-sided sampled surface deviation (approximate hausdorff/mean)."""
    step_mesh = step_to_mesh(step_path)
    if len(step_mesh.faces) == 0 or len(mesh.faces) == 0:
        raise ValueError("empty geometry — nothing to compare")
    a = mesh.sample(samples)
    b = step_mesh.sample(samples)
    d_ab = trimesh.proximity.closest_point(step_mesh, a)[1]
    d_ba = trimesh.proximity.closest_point(mesh, b)[1]
    return {
        "max": float(max(d_ab.max(), d_ba.max())),
        "mean": float((d_ab.mean() + d_ba.mean()) / 2),
    }

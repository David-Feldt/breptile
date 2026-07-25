"""Load and repair input meshes."""

from __future__ import annotations

import numpy as np
import trimesh


class MeshError(RuntimeError):
    pass


def load_mesh(path: str, force: bool = False, max_triangles: int | None = None) -> trimesh.Trimesh:
    """Load an STL (or any trimesh-supported mesh), repair it, and validate.

    Raises MeshError if the mesh is not watertight and ``force`` is False.
    """
    mesh = trimesh.load(path, force="mesh", process=True)
    if not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) == 0:
        raise MeshError(f"{path}: no triangle geometry found")

    mesh.remove_unreferenced_vertices()
    mesh.update_faces(mesh.nondegenerate_faces())
    trimesh.repair.fix_normals(mesh)

    if not mesh.is_watertight:
        trimesh.repair.fill_holes(mesh)

    if not mesh.is_watertight:
        # Last resort: rebuild through manifold3d, which can heal many broken meshes.
        try:
            import manifold3d

            m = manifold3d.Manifold(
                manifold3d.Mesh(
                    vert_properties=np.asarray(mesh.vertices, dtype=np.float32),
                    tri_verts=np.asarray(mesh.faces, dtype=np.uint32),
                )
            )
            out = m.to_mesh()
            healed = trimesh.Trimesh(
                vertices=np.asarray(out.vert_properties)[:, :3],
                faces=np.asarray(out.tri_verts),
                process=True,
            )
            if healed.is_watertight and len(healed.faces) > 0:
                mesh = healed
        except Exception:
            pass

    if not mesh.is_watertight and not force:
        raise MeshError(
            f"{path}: mesh is not watertight after repair "
            f"({len(mesh.faces)} faces). Re-run with --force to convert as an open shell."
        )

    if max_triangles and len(mesh.faces) > max_triangles:
        mesh = mesh.simplify_quadric_decimation(face_count=max_triangles)
        trimesh.repair.fix_normals(mesh)

    return mesh

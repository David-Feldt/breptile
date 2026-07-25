"""Render a before/after grid: STL wireframe (left half of each cell) vs
converted STEP with BREP edges (right half)."""
import json
import os

import numpy as np
import trimesh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

from stl2step.verify import step_to_mesh

HERE = os.path.dirname(os.path.abspath(__file__))


def step_edges(step_path, deflection):
    from build123d import import_step
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_ShapeEnum
    from OCP.TopoDS import TopoDS
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.GCPnts import GCPnts_QuasiUniformDeflection

    shape = import_step(step_path).wrapped
    lines = []
    exp = TopExp_Explorer(shape, TopAbs_ShapeEnum.TopAbs_EDGE)
    while exp.More():
        try:
            curve = BRepAdaptor_Curve(TopoDS.Edge_s(exp.Current()))
            disc = GCPnts_QuasiUniformDeflection(curve, deflection)
            if disc.IsDone():
                pts = np.array([[disc.Value(i).X(), disc.Value(i).Y(), disc.Value(i).Z()]
                                for i in range(1, disc.NbPoints() + 1)])
                if len(pts) >= 2:
                    lines.append(pts)
        except Exception:
            pass
        exp.Next()
    return lines


def draw(ax, mesh, title, edge_lines=None, facecolor="#dfe7ef", edgecolor="#5b7ea6"):
    tris = mesh.vertices[mesh.faces]
    if edge_lines is None:
        ax.add_collection3d(Poly3DCollection(tris, facecolor=facecolor,
                                             edgecolor=edgecolor, linewidth=0.2))
    else:
        ax.add_collection3d(Poly3DCollection(tris, facecolor="#ece5d5", edgecolor="none"))
        ax.add_collection3d(Line3DCollection(edge_lines, colors="#7a5c22", linewidths=0.7))
    c, e = mesh.bounding_box.centroid, mesh.extents.max() * 0.6
    ax.set_xlim(c[0]-e, c[0]+e); ax.set_ylim(c[1]-e, c[1]+e); ax.set_zlim(c[2]-e, c[2]+e)
    ax.set_axis_off()
    ax.view_init(elev=35, azim=-60)
    ax.set_title(title, fontsize=8)


results = {r["model"]: r for r in json.load(open(os.path.join(HERE, "results.json")))}
models = [m for m, r in results.items() if r.get("status") == "ok" and r.get("valid")]
models.sort()

cols = 6
rows = int(np.ceil(len(models) * 2 / cols))
fig = plt.figure(figsize=(3.0 * cols, 3.0 * rows), dpi=110)

i = 0
for name in models:
    r = results[name]
    stl = os.path.join(HERE, "corpus", name)
    step = os.path.join(HERE, "out", os.path.splitext(name)[0] + ".step")
    try:
        mesh = trimesh.load(stl, force="mesh")
        defl = mesh.extents.max() / 400
        smesh = step_to_mesh(step, defl)
        edges = step_edges(step, defl)
    except Exception as e:
        print("skip", name, e)
        continue
    ax1 = fig.add_subplot(rows, cols, i + 1, projection="3d")
    draw(ax1, mesh, f"{name}\nSTL · {r['triangles']:,} tris")
    ax2 = fig.add_subplot(rows, cols, i + 2, projection="3d")
    label = f"STEP · {r['faces']} faces · {r.get('analytic_pct', 0)}% analytic"
    draw(ax2, smesh, label, edge_lines=edges)
    i += 2
    print("done", name, flush=True)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "benchmark_grid.png"), bbox_inches="tight")
print("saved grid")

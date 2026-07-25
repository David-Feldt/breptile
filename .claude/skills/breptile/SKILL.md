---
name: breptile
description: Convert an STL mesh to clean editable STEP using the breptile pipeline, then inspect the report and manually rebuild regions the automatic fitter could not handle (hybrid mesh-to-BREP workflow).
---

# breptile hybrid conversion

Convert STL → STEP with the `breptile` CLI, then act as the second stage: fix what the
algorithm couldn't. In this repo, run everything with `.venv/bin/python` /
`.venv/bin/breptile`. Elsewhere, `pip install breptile` and use `breptile` / `python`
directly.

## Step 1 — automatic pass

```bash
.venv/bin/breptile input.stl output.step --report report.json --verify
```

Read `report.json`:
- `segmentation` / `regions`: counts of plane/cylinder/sphere/freeform regions and how many
  fell back to tessellation (`fallback`, `freeform_triangles`).
- `deviation.max`: sampled two-sided deviation between the STEP and the input mesh. Compare
  against `tolerance` — values at the input mesh's chord error are expected and fine.
- `valid_brep`: must be true.

If `freeform_triangles == 0`, `fallback == 0`, and the deviation is acceptable: done.

## Step 2 — diagnose failures

Common causes and knobs:
- **Coarse tessellation** (facet angles > 30°): re-run with a looser smooth angle — call the
  Python API: `convert(inp, out, mode="fit", tol=...)`, or pre-inspect with
  `breptile.segment.segment(mesh, tol, smooth_angle=np.radians(45))`.
- **Tolerance too tight for a noisy mesh**: raise `--tol` (fits are rejected when max
  deviation of region vertices exceeds it).
- **Partial cylinders / fillets / spheres**: v1 only rebuilds full-wrap cylinders with two
  circular rims; everything else tessellates. These are LLM-rebuild candidates.

To see what a region looks like, render it:

```python
import trimesh, numpy as np
from breptile.mesh import load_mesh
from breptile.segment import segment
mesh = load_mesh("input.stl")
regions = segment(mesh, tol)
bad = [r for r in regions if r.kind == "freeform"]
sub = mesh.submesh([bad[0].faces], append=True)
sub.export("region0.glb")  # or scene screenshot via trimesh
```

## Step 3 — LLM rebuild of failed regions

When the part's design intent is recognizable, rebuild it as build123d code instead of
patching facets:

1. Measure from the mesh: bounding box, fitted primitive params (`Region.params` has
   axis/center/radius even for rejected fits), cross-sections
   (`mesh.section(plane_origin=..., plane_normal=...)`).
2. Write a build123d script reproducing the part; export with `export_step(part, "output.step")`.
3. **Never snap dimensions** unless the user asks — this project's contract is tight
   tolerance to the mesh. Use measured values.

Alternatively, boolean-patch: import the auto-converted STEP with `import_step`, cut away
the bad region's bounding volume, and fuse a cleanly modeled replacement.

## Step 4 — verify

Always finish by measuring the result against the original mesh:

```python
from breptile.verify import deviation
from breptile.mesh import load_mesh
print(deviation("output.step", load_mesh("input.stl")))
```

`max` must be ≤ the agreed tolerance (or the input's chord error). Also confirm
`valid_brep` via `breptile.brep.is_valid` and, for editability, that expected features are
analytic (see `tests/test_convert.py::_cyl_face_count` for the pattern).

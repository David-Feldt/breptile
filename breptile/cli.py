"""breptile command-line interface."""

from __future__ import annotations

import argparse
import json
import sys

from breptile.mesh import MeshError
from breptile.pipeline import convert


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="breptile",
        description="Convert an STL mesh to a clean, editable STEP BREP model.",
    )
    p.add_argument("input", help="input mesh file (STL/OBJ/3MF/PLY)")
    p.add_argument("output", help="output STEP file")
    p.add_argument(
        "--mode",
        choices=["auto", "tessellated", "prismatic", "fit"],
        default="auto",
        help="auto: fit primitives with per-region fallback (default); "
        "prismatic: planes only; tessellated: one face per triangle (always succeeds)",
    )
    p.add_argument("--tol", type=float, default=None,
                   help="max surface deviation (default: bbox diagonal * 1e-4)")
    p.add_argument("--schema", choices=["AP214", "AP242"], default="AP214")
    p.add_argument("--report", metavar="JSON", help="write conversion report to this file")
    p.add_argument("--max-triangles", type=int, default=None,
                   help="decimate input above this triangle count")
    p.add_argument("--force", action="store_true",
                   help="convert non-watertight meshes as open shells")
    p.add_argument("--verify", action="store_true",
                   help="re-tessellate the STEP and report deviation vs input")
    args = p.parse_args(argv)

    try:
        report = convert(
            args.input,
            args.output,
            mode="fit" if args.mode == "auto" else args.mode,
            tol=args.tol,
            schema=args.schema,
            force=args.force,
            max_triangles=args.max_triangles,
        )
    except MeshError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.verify:
        from breptile.mesh import load_mesh
        from breptile.verify import deviation

        try:
            mesh = load_mesh(args.input, force=True)
            report["deviation"] = deviation(args.output, mesh)
        except Exception as e:
            report["deviation"] = {"error": f"verification failed: {e}"}

    if args.report:
        with open(args.report, "w") as f:
            json.dump(report, f, indent=2, default=str)

    print(json.dumps(report, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())

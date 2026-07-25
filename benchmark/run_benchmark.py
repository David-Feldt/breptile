"""Benchmark stl2step across a corpus of models. Each conversion runs in a
subprocess with a timeout so one pathological mesh can't hang the sweep."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import urllib.request

BASE = "https://github.com/mikedh/trimesh/raw/main/models/"
MODELS = [
    "20mm-xyz-cube.stl",
    "7_8ths_cube.stl",
    "angle_block.STL",
    "box.STL",
    "busted.STL",
    "cylinder.stl",
    "featuretype.STL",
    "idler_riser.STL",
    "octagonal_pocket.stl",
    "origin_inside.STL",
    "plate_holes.STL",
    "round.stl",
    "soup.stl",
    "teapot.stl",
    "torus.STL",
    "unit_sphere.STL",
    "ADIS16480.STL",
    "1002_tray_bottom.STL",
]

HERE = os.path.dirname(os.path.abspath(__file__))
CORPUS = os.path.join(HERE, "corpus")
OUT = os.path.join(HERE, "out")
TIMEOUT = 240


def fetch(name: str) -> str:
    os.makedirs(CORPUS, exist_ok=True)
    path = os.path.join(CORPUS, name)
    if not os.path.exists(path):
        urllib.request.urlretrieve(BASE + name, path)
    return path


def run_one(name: str) -> dict:
    stl = fetch(name)
    step = os.path.join(OUT, os.path.splitext(name)[0] + ".step")
    report = step + ".json"
    row: dict = {"model": name, "stl_kb": round(os.path.getsize(stl) / 1024)}

    def attempt(extra_args):
        cmd = [sys.executable, "-m", "stl2step.cli", stl, step,
               "--report", report, "--verify", *extra_args]
        t0 = time.time()
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        return p, time.time() - t0

    try:
        p, dt = attempt([])
        if p.returncode == 2 and "not watertight" in p.stderr:
            row["note"] = "not watertight → --force"
            p, dt = attempt(["--force"])
        if p.returncode != 0:
            row["status"] = "error"
            row["error"] = (p.stderr or p.stdout).strip().splitlines()[-1][:200]
            return row
        r = json.load(open(report))
        seg = r.get("segmentation", {})
        reg = r.get("regions", {})
        row.update(
            status="ok",
            seconds=round(dt, 1),
            triangles=r["triangles"],
            faces=r["faces"],
            planes=seg.get("plane", 0),
            cylinders=seg.get("cylinder", 0),
            spheres=seg.get("sphere", 0),
            freeform=seg.get("freeform", 0),
            fallback=reg.get("fallback", 0),
            freeform_tris=reg.get("freeform_triangles", 0),
            valid=r["valid_brep"],
            dev_max=round(r["deviation"]["max"], 5) if "max" in r.get("deviation", {}) else None,
            dev_mean=round(r["deviation"]["mean"], 6) if "mean" in r.get("deviation", {}) else None,
            step_kb=round(os.path.getsize(step) / 1024),
        )
        # Analytic ratio: how much of the mesh got promoted to analytic faces.
        if r["triangles"]:
            row["analytic_pct"] = round(100 * (1 - row["freeform_tris"] / r["triangles"]))
    except subprocess.TimeoutExpired:
        row["status"] = "timeout"
        row["seconds"] = TIMEOUT
    except Exception as e:
        row["status"] = "error"
        row["error"] = str(e)[:200]
    return row


def main():
    os.makedirs(OUT, exist_ok=True)
    results = []
    for name in MODELS:
        print(f"=== {name}", flush=True)
        row = run_one(name)
        print(json.dumps(row), flush=True)
        results.append(row)
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("wrote results.json")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate a time-dependent FlatLand case: orbiting view + evolving scalar field.

Dependency-free: parses the OBJ itself and uses only the Python standard library
(no numpy / pyvista), so the example reproduces with a stock Python 3.

It writes a single multi-column field file (one row per vertex, one column per
timestep) and a batch file whose lines select a column per timestep:

    a <azimuth> <elevation> <resolution> field.txt@<t>

The whole time series therefore lives in two files regardless of step count.

Usage:
    python3 generate_field.py <mesh.obj> <out_dir> [num_steps] [resolution]

Example:
    python3 generate_field.py bunny.obj data 1000 0.002
    ../../src/build/flatland bunny.obj -b data/timeseries.txt -j > results.json
"""
import math
import os
import sys


def load_vertices(path):
    verts = []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                _, x, y, z = line.split()[:4]
                verts.append((float(x), float(y), float(z)))
    if not verts:
        raise SystemExit(f"no vertices found in {path}")
    return verts


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    mesh = sys.argv[1]
    out_dir = sys.argv[2]
    steps = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    res = float(sys.argv[4]) if len(sys.argv) > 4 else 0.001

    verts = load_vertices(mesh)
    os.makedirs(out_dir, exist_ok=True)

    # Normalize x into [0, 1] so the traveling wave is mesh-independent.
    xs = [v[0] for v in verts]
    x0, x1 = min(xs), max(xs)
    span = (x1 - x0) or 1.0
    nx = [(v[0] - x0) / span for v in verts]

    # One field column per timestep; a wave that sweeps along x as time advances.
    freq = 4.0
    field_path = os.path.join(out_dir, "field.txt")
    with open(field_path, "w") as ff:
        for i in range(len(verts)):
            row = (f"{math.sin(freq * math.pi * nx[i] - 2.0 * math.pi * t / steps):.6f}"
                   for t in range(steps))
            ff.write(" ".join(row) + "\n")

    # Batch: orbit the mesh in azimuth (slight tilt), one column per timestep.
    # Field path is a bare name; FlatLand resolves it against the batch file's dir.
    batch_path = os.path.join(out_dir, "timeseries.txt")
    with open(batch_path, "w") as batch:
        batch.write("# a <azimuth> <elevation> <resolution> field.txt@<timestep>\n")
        for t in range(steps):
            az = 360.0 * t / steps
            batch.write(f"a {az:.4f} 20 {res} field.txt@{t}\n")

    print(f"Wrote {steps}-column field to {field_path} and batch to {batch_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate a time-dependent FlatLand case: rotating view + evolving scalar field.

Dependency-free: parses the OBJ itself and uses only the Python standard library
(no numpy / pyvista). This makes the example reproducible with a stock Python 3.

For each timestep it writes one node-field file (one value per vertex) and appends
a line to a batch file of the form expected by `flatland -b`:

    <nx> <ny> <nz> <resolution> <field_file>

Usage:
    python3 generate_field.py <mesh.obj> <out_dir> [num_steps] [resolution]

Example:
    python3 generate_field.py bunny.obj data 200 0.001
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

    batch_path = os.path.join(out_dir, "timeseries.txt")
    with open(batch_path, "w") as batch:
        batch.write("# nx ny nz resolution field_file  (one timestep per line)\n")
        for t in range(steps):
            # Field: a wave that sweeps along x as time advances.
            phase = 2.0 * math.pi * t / steps
            freq = 4.0
            # Field file lives beside the batch file; FlatLand resolves the batch's
            # relative data paths against the batch file's directory, so the batch
            # only needs the bare filename and the case stays portable.
            field_name = f"field_{t:04d}.txt"
            with open(os.path.join(out_dir, field_name), "w") as ff:
                ff.write("\n".join(
                    f"{math.sin(freq * math.pi * ((v[0] - x0) / span) - phase):.6f}"
                    for v in verts
                ))
                ff.write("\n")

            # View: orbit the mesh about the z-axis with a slight tilt.
            ang = 2.0 * math.pi * t / steps
            nx, ny, nz = math.cos(ang), math.sin(ang), 0.3
            batch.write(f"{nx:.6f} {ny:.6f} {nz:.6f} {res} {field_name}\n")

    print(f"Wrote {steps} timesteps to {batch_path}")


if __name__ == "__main__":
    main()

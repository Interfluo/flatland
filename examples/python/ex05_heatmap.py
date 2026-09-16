#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
Example 5 — Getting at the raster itself.

The question this answers: *I do not just want the number, I want to see where it
came from.* `render()` hands back the coverage mask and the interpolated
per-pixel field values, so you can plot them, post-process them, or feed them
into something else.

    python3 examples/python/ex05_heatmap.py [outdir]

Covered:
  - the coverage mask, and checking it against covered_pixels
  - per-pixel field values, and recomputing the integral from them by hand
  - writing a PPM, and converting to PNG with no image library
  - matplotlib, if it happens to be installed
  - row order: row 0 is the BOTTOM row in mesh space
"""

import math
import os
import struct
import sys
import zlib

from _setup import load, mesh_path, rule

fl = load()


def write_png(path, width, height, rgb_rows):
    """Minimal PNG writer — zlib is in the standard library, so no dependency.

    rgb_rows is a list of `height` bytearrays, each 3*width bytes, TOP row first.
    """
    raw = b"".join(b"\x00" + bytes(row) for row in rgb_rows)

    def chunk(tag, data):
        c = struct.pack(">I", len(data)) + tag + data
        return c + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)

    header = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    with open(path, "wb") as fh:
        fh.write(b"\x89PNG\r\n\x1a\n")
        fh.write(chunk(b"IHDR", header))
        fh.write(chunk(b"IDAT", zlib.compress(raw, 9)))
        fh.write(chunk(b"IEND", b""))


def ramp(t):
    """Blue -> cyan -> yellow -> orange -> red, matching FlatLand's own PPM."""
    stops = [(0, 0, 1), (0, 1, 1), (1, 1, 0), (1, .5, 0), (1, 0, 0)]
    t = max(0.0, min(1.0, t)) * 4.0
    i = min(3, int(t))
    f = t - i
    return tuple(int(((1 - f) * stops[i][c] + f * stops[i + 1][c]) * 255)
                 for c in range(3))


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(outdir, exist_ok=True)
    print("FlatLand %s   writing into %s" % (fl.__version__, os.path.abspath(outdir)))

    mesh = fl.Mesh.load(mesh_path("sphere_areas", "sphere_fine.obj"))
    view = (0.4, 0.3, 1.0)

    # f = cos(incidence): brightest facing the camera, falling to zero at the limb.
    verts = mesh.vertices
    n = math.sqrt(sum(c * c for c in view))
    unit_view = tuple(c / n for c in view)
    field = [-sum(v[k] * unit_view[k] for k in range(3)) /
             math.sqrt(sum(c * c for c in v)) for v in verts]

    # ---------------------------------------------------------------------
    rule("Rendering")

    img = mesh.render(view, field=field, resolution=4e-3, precision="double")
    r = img.result
    print("  raster      %d x %d" % (img.width, img.height))
    print("  covered     %d pixels of %d" % (r.covered_pixels, img.width * img.height))
    print("  area        %.6f" % r.area)
    print("  field range %.4f .. %.4f, mean %.4f" % (r.min, r.max, r.average))

    # ---------------------------------------------------------------------
    rule("The mask")

    # mask is width*height bytes, non-zero where a pixel is covered. Index it
    # as [y * width + x]. Row 0 is the BOTTOM row in mesh space, which is the
    # opposite of most image formats — flip when displaying.
    mask = img.mask
    covered = sum(1 for m in mask if m)
    print("  mask length %d  (= width x height: %s)"
          % (len(mask), len(mask) == img.width * img.height))
    print("  covered in the mask: %d  (matches covered_pixels: %s)"
          % (covered, covered == r.covered_pixels))

    # ---------------------------------------------------------------------
    rule("The values, and recomputing the integral by hand")

    # values holds the interpolated field at each covered pixel. Summing it and
    # multiplying by the pixel area reproduces the integral FlatLand reported —
    # which is exactly how it is defined.
    values = img.values
    pixel_area = 4e-3 * 4e-3
    manual = sum(v for v, m in zip(values, mask) if m) * pixel_area
    print("  FlatLand's integral : %.8f" % r.integral)
    print("  recomputed by hand  : %.8f" % manual)
    print("  agree to 1e-9       : %s" % (abs(manual - r.integral) < 1e-9))

    # ---------------------------------------------------------------------
    rule("Writing images")

    ppm = os.path.join(outdir, "sphere_heatmap.ppm")
    img.save_ppm(ppm)
    print("  PPM  %s  (%d bytes)" % (ppm, os.path.getsize(ppm)))

    # The same raster as a PNG, built here so the example needs no image library.
    lo, hi = r.min, r.max
    span = (hi - lo) or 1.0
    rows = []
    for y in range(img.height - 1, -1, -1):          # flip: PNG is top-down
        row = bytearray()
        for x in range(img.width):
            i = y * img.width + x
            if mask[i]:
                row += bytes(ramp((values[i] - lo) / span))
            else:
                row += bytes((30, 30, 35))
        rows.append(row)
    png = os.path.join(outdir, "sphere_heatmap.png")
    write_png(png, img.width, img.height, rows)
    print("  PNG  %s  (%d bytes)" % (png, os.path.getsize(png)))

    # ---------------------------------------------------------------------
    rule("matplotlib")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib is not installed here. With it you would write:")
        print("    import numpy as np, matplotlib.pyplot as plt")
        print("    v = np.array(img.values).reshape(img.height, img.width)")
        print("    v[~np.array(img.mask, bool).reshape(v.shape)] = np.nan")
        print("    plt.imshow(v, origin='lower'); plt.colorbar(); plt.show()")
    else:
        grid = [[values[y * img.width + x] if mask[y * img.width + x] else float("nan")
                 for x in range(img.width)] for y in range(img.height)]
        fig, ax = plt.subplots(figsize=(5, 5))
        im = ax.imshow(grid, origin="lower")          # origin='lower' for row 0 = bottom
        ax.set_axis_off()
        fig.colorbar(im, ax=ax, label="cos(incidence)")
        out = os.path.join(outdir, "sphere_matplotlib.png")
        fig.savefig(out, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print("  wrote %s" % out)

    # ---------------------------------------------------------------------
    rule("An empty view")

    # A view that covers nothing yields a valid 0x0 image rather than a stale
    # raster from a previous call, and refuses to write a malformed file.
    flat_v = [(0, 0, 0), (1, 0, 0), (1, 1, 0)]
    flat_f = [(0, 1, 2)]
    with fl.Mesh(flat_v, flat_f) as flat:
        edge_on = flat.render((1, 0, 0), resolution=1e-2)
        print("  edge-on raster: %dx%d, mask length %d"
              % (edge_on.width, edge_on.height, len(edge_on.mask)))
        try:
            edge_on.save_ppm(os.path.join(outdir, "should_not_exist.ppm"))
            print("  ERROR: it should have refused to write that")
        except fl.FlatlandError as exc:
            print("  save_ppm refused, as it should: %s" % exc)

    img.close()
    mesh.close()


if __name__ == "__main__":
    main()

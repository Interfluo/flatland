#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""
Example 5 — Getting at the raster itself.

The question this answers: *I do not just want the number, I want to see where it
came from, or feed it into something else.* `render()` hands back the coverage
mask and the interpolated per-pixel field values, and the raster can be written
two ways:

    save_png()   a picture. The field is ramped to 8 bits per channel, which is
                 right for looking at and wrong for computing with.
    save_npy()   the numbers. float64, NaN off the silhouette, np.load()-able.

Both are produced by the C library itself — FlatLand links no image or
compression library, so the PNG encoder is built in.

    python3 examples/python/ex05_heatmap.py [outdir]

Covered:
  - the coverage mask, and checking it against covered_pixels
  - recomputing the integral from the raster by hand
  - PNG for looking at, NPY for computing with
  - matplotlib, if it happens to be installed
  - row order: row 0 is the BOTTOM row in mesh space
"""

import math
import os
import sys

from _setup import load, mesh_path, rule

fl = load()


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(outdir, exist_ok=True)
    print("FlatLand %s   writing into %s" % (fl.__version__, os.path.abspath(outdir)))

    mesh = fl.Mesh.load(mesh_path("sphere_areas", "sphere_fine.obj"))
    view = (0.4, 0.3, 1.0)
    res = 4e-3

    # f = cos(incidence): brightest facing the camera, falling to zero at the limb.
    verts = mesh.vertices
    n = math.sqrt(sum(c * c for c in view))
    unit_view = tuple(c / n for c in view)
    field = [-sum(v[k] * unit_view[k] for k in range(3)) /
             math.sqrt(sum(c * c for c in v)) for v in verts]

    # ---------------------------------------------------------------------
    rule("Rendering")

    img = mesh.render(view, field=field, resolution=res, precision="double")
    r = img.result
    print("  raster      %d x %d" % (img.width, img.height))
    print("  covered     %d pixels of %d" % (r.covered_pixels, img.width * img.height))
    print("  area        %.6f" % r.area)
    print("  field range %.4f .. %.4f, mean %.4f" % (r.min, r.max, r.average))

    # ---------------------------------------------------------------------
    rule("The mask")

    # mask is width*height bytes, non-zero where a pixel is covered. Index it as
    # [y * width + x]. Row 0 is the BOTTOM row in mesh space, the opposite of
    # most image formats — save_png flips on the way out, save_npy does not.
    mask = img.mask
    covered = sum(1 for m in mask if m)
    print("  mask length %d  (= width x height: %s)"
          % (len(mask), len(mask) == img.width * img.height))
    print("  covered in the mask: %d  (matches covered_pixels: %s)"
          % (covered, covered == r.covered_pixels))

    # ---------------------------------------------------------------------
    rule("Recomputing the integral by hand")

    # values holds the interpolated field at each covered pixel. Summing it and
    # multiplying by the pixel area reproduces the integral FlatLand reported —
    # which is exactly how the integral is defined.
    values = img.values
    pixel_area = res * res
    manual = sum(v for v, m in zip(values, mask) if m) * pixel_area
    print("  FlatLand's integral : %.8f" % r.integral)
    print("  recomputed by hand  : %.8f" % manual)
    print("  agree to 1e-9       : %s" % (abs(manual - r.integral) < 1e-9))

    # ---------------------------------------------------------------------
    rule("PNG — for looking at")

    png = os.path.join(outdir, "sphere_heatmap.png")
    img.save_png(png)
    raw_rgb = img.width * img.height * 3
    size = os.path.getsize(png)
    print("  %s  (%d bytes)" % (png, size))
    print("  %.1fx smaller than the same raster uncompressed (%d bytes)"
          % (raw_rgb / size, raw_rgb))
    print("  the colour ramp quantises to 8 bits per channel, so this is a")
    print("  picture of the field, not the field")

    # ---------------------------------------------------------------------
    rule("NPY — for computing with")

    npy = os.path.join(outdir, "sphere_field.npy")
    img.save_npy(npy)
    print("  %s  (%d bytes, float64)" % (npy, os.path.getsize(npy)))
    print("  shape (%d, %d), NaN wherever nothing was covered" % (img.height, img.width))
    print("  row 0 is the BOTTOM row, so plot it with origin='lower'")

    if fl.HAS_NUMPY:
        import numpy as np
        arr = np.load(npy)
        cov = ~np.isnan(arr)
        print("  reloaded: shape %s, %d covered, mean %.10f"
              % (arr.shape, cov.sum(), arr[cov].mean()))
        print("  matches the reported mean: %s"
              % (abs(arr[cov].mean() - r.average) < 1e-9))
        print("  integral from the array  : %.8f"
              % (arr[cov].sum() * pixel_area))
    else:
        print("  (NumPy is not installed here, so it is not reloaded — but the")
        print("   file is written by the C library and needs no NumPy to produce)")

    # ---------------------------------------------------------------------
    rule("matplotlib")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib is not installed here. With it, the .npy makes this")
        print("  a two-liner:")
        print("    a = np.load('sphere_field.npy')            # NaN off the silhouette")
        print("    plt.imshow(a, origin='lower'); plt.colorbar()")
    else:
        if fl.HAS_NUMPY:
            import numpy as np
            grid = np.load(npy)
        else:
            grid = [[values[y * img.width + x] if mask[y * img.width + x] else float("nan")
                     for x in range(img.width)] for y in range(img.height)]
        fig, ax = plt.subplots(figsize=(5, 5))
        im = ax.imshow(grid, origin="lower")          # row 0 is the bottom row
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
        for name, fn in (("save_png", edge_on.save_png), ("save_npy", edge_on.save_npy)):
            try:
                fn(os.path.join(outdir, "should_not_exist"))
                print("  ERROR: %s should have refused" % name)
            except fl.FlatlandError as exc:
                print("  %s refused, as it should: %s" % (name, exc))

    img.close()
    mesh.close()


if __name__ == "__main__":
    main()

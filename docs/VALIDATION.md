# Validation

FlatLand is checked against closed-form results, not against its own previous
output. Every expected value on this page is derived analytically; nothing is a
stored baseline, so a wrong answer cannot be blessed into a regression fixture.

Reproduce everything here with:

```shell
make validate            # the full study, ~3 s
make test                # the quick subset, as part of the test suite
```

The numbers below were produced on Linux/x86-64 with GCC 13 in double precision.
They will shift in the last digit or two on other toolchains; the tolerances are
set well clear of that.

---

## Method

Two error sources have to be kept apart, because they behave completely
differently and conflating them hides both:

- **Rasterization error** — FlatLand decides coverage by testing whether a pixel
  *centre* falls inside a triangle. A silhouette therefore has a boundary band
  roughly one pixel wide that is either counted or not, giving an error that
  shrinks with the pixel size `h`.
- **Mesh error** — a triangle mesh of a curved body is not the body. An
  icosphere is an *inscribed* polyhedron, so it genuinely projects to less than
  `πr²`. That is a property of the input, not a defect in the tool.

The trick that separates them is that a convex polyhedron's projected area has a
closed form. So for a faceted sphere we can compute what the tool *should*
report for that exact mesh, independently of what a true sphere would give.

### The three results everything rests on

**1. Projected area of a convex polyhedron.** Every line through a convex body
crosses its surface exactly twice, so each direction sees precisely half the
total projected face area:

$$A(\hat n) = \tfrac{1}{2}\sum_i A_i\,\lvert \hat n \cdot \hat n_i\rvert$$

With $A_i = \tfrac12\lVert c_i\rVert$ and $\hat n_i = c_i/\lVert c_i\rVert$ for
the triangle cross product $c_i = (v_1-v_0)\times(v_2-v_0)$, this collapses to a
form with no normalization at all:

$$A(\hat n) = \tfrac{1}{4}\sum_i \lvert \hat n \cdot c_i \rvert$$

**2. Cauchy's projection formula.** Averaged over directions uniform on the
sphere, a convex body's projected area is a quarter of its surface area:

$$\langle A\rangle = S/4$$

This is an integral identity, so it probes many directions at once rather than a
few hand-picked ones.

**3. Mean of a linear field over a triangle.** Barycentric interpolation
reproduces a linear function exactly, and the mean of a linear function over a
triangle is the mean of its vertex values:

$$\frac{1}{A}\int f\,dA = \frac{f_0+f_1+f_2}{3}$$

The formulas themselves are checked before anything else runs
(`validation/self_check.py`, 33 identities) against values known from elsewhere —
that a cube down its body diagonal projects to a regular hexagon of area
$\sqrt3$, that a regular tetrahedron has $S=\sqrt3 e^2$, that an inscribed
polyhedron must understate the sphere and improve monotonically with
subdivision. Without that step a mistake in the expectations would make every
study agree with the wrong answer and report success.

---

## Results

### Exact convex polyhedra — isolating rasterization error

These meshes *are* the shape, so the entire discrepancy is rasterization.
Pixel size `h = 5×10⁻⁴`.

| Shape | View | FlatLand | Exact | Rel. error |
| :--- | :--- | ---: | ---: | ---: |
| unit cube | (1,0,0) | 0.99999714 | 1.00000000 | 2.9×10⁻⁶ |
| unit cube | (1,1,0) | 1.41400015 | 1.41421356 | 1.5×10⁻⁴ |
| unit cube | (1,1,1) | 1.73187971 | 1.73205081 | 9.9×10⁻⁵ |
| unit cube | (0.3,0.9,−0.31) | 1.51280737 | 1.51295314 | 9.6×10⁻⁵ |
| tetrahedron | (1,0,0) | 0.49984181 | 0.50000000 | 3.2×10⁻⁴ |
| tetrahedron | (1,1,1) | 0.43301481 | 0.43301270 | 4.9×10⁻⁶ |
| octahedron | (0,0,1) | 1.99994552 | 2.00000000 | 2.7×10⁻⁵ |
| octahedron | (1,1,1) | 1.73187149 | 1.73205081 | 1.0×10⁻⁴ |

21 direction/shape combinations, worst relative error **3.2×10⁻⁴**.

The unit cube also has the especially simple closed form
$A(\hat n) = \lvert n_x\rvert + \lvert n_y\rvert + \lvert n_z\rvert$. Over 40
quasi-uniform directions the worst relative error is **3.5×10⁻⁴**.

### Cauchy's identity — many directions at once

400 quasi-uniform directions, `h = 2×10⁻³`.

| Shape | Mean projected area | S/4 | Rel. error |
| :--- | ---: | ---: | ---: |
| unit cube | 1.50001838 | 1.50000000 | 1.2×10⁻⁵ |
| tetrahedron | 0.43301730 | 0.43301270 | 1.1×10⁻⁵ |
| octahedron | 1.73205246 | 1.73205081 | 9.6×10⁻⁷ |
| icosphere(3) | 3.12663817 | 3.12662318 | 4.8×10⁻⁶ |

Averaging cancels the per-direction boundary error almost completely, which is
why these are an order of magnitude tighter than the single-direction table.

### Sphere — the two error sources, separated

Radius 1, viewed along the body diagonal, `h = 10⁻³`.

| Subdiv | Triangles | FlatLand | Exact polyhedron | Raster error | vs. true sphere |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 80 | 2.9467282 | 2.9467408 | 4.3×10⁻⁶ | 6.2×10⁻² |
| 2 | 320 | 3.0849884 | 3.0849753 | 4.2×10⁻⁶ | 1.8×10⁻² |
| 3 | 1280 | 3.1281223 | 3.1281120 | 3.3×10⁻⁶ | 4.3×10⁻³ |
| 4 | 5120 | 3.1380982 | 3.1380956 | 8.5×10⁻⁷ | 1.1×10⁻³ |
| 5 | 20480 | 3.1407144 | 3.1407082 | 2.0×10⁻⁶ | 2.8×10⁻⁴ |

This is the most informative table on the page. The **raster error column is flat
at ~10⁻⁶** and does not improve with subdivision, because it has nothing to do
with the mesh. The **final column falls by a factor of four per subdivision**,
which is the inscribed polyhedron converging on the sphere as the square of the
edge length. FlatLand reproduces whatever mesh it is handed to about six digits;
the remaining gap to `π` is the mesh's, and the only way to close it is a finer
mesh, not a finer raster.

### Cylinder — curved and flat surfaces together

$A(\theta) = 2rh\sin\theta + \pi r^2\cos\theta$, with $r=1$, $h=2$, a 512-gon
prism, `h = 10⁻³`.

| θ | FlatLand | Analytic | Rel. error |
| ---: | ---: | ---: | ---: |
| 0° | 3.14144731 | 3.14159265 | 4.6×10⁻⁵ |
| 30° | 4.72061968 | 4.72069905 | 1.7×10⁻⁵ |
| 45° | 5.04979467 | 5.04986859 | 1.5×10⁻⁵ |
| 60° | 5.03486061 | 5.03489794 | 7.4×10⁻⁶ |
| 90° | 4.00000048 | 4.00000000 | 1.2×10⁻⁷ |

Both terms contribute at intermediate angles, so this exercises the end caps and
the curved wall simultaneously.

### Linear fields — exact interpolation

`h = 5×10⁻⁴`, double precision.

| Case | Quantity | FlatLand | Exact | Rel. error |
| :--- | :--- | ---: | ---: | ---: |
| unit square, f = x | ∫f dA | 0.49999867 | 1/2 | 2.7×10⁻⁶ |
| unit square, f = x | mean | 0.50001529 | 1/2 | 3.1×10⁻⁵ |
| right triangle, f = x | ∫f dA | 0.16660307 | 1/6 | 3.8×10⁻⁴ |
| right triangle, f = x | mean | 0.33320763 | 1/3 | 3.8×10⁻⁴ |
| unit cube, f = y along +x | ∫f dA | 0.49999867 | 1/2 | 2.7×10⁻⁶ |

The cube case is doubly useful: viewed along +x the visible face is `x = 0`, so
it only produces 1/2 if the depth test resolves to the **near** surface. An
inverted backface cull gives the same *area* and a different *mean*, which is
exactly the failure mode that went undetected before this suite existed.

### Lambertian sphere — a non-trivial field with a closed form

$f = \cos\theta$ between the surface normal and the view direction. At projected
radius $\rho$ the cosine is $\sqrt{1-(\rho/r)^2}$, so

$$\int f\,dA = \int_0^r \sqrt{1-(\rho/r)^2}\;2\pi\rho\,d\rho = \tfrac{2}{3}\pi r^2,
\qquad \langle f\rangle = \tfrac{2}{3}$$

| Subdiv | Triangles | ∫f dA | Exact | Mean | Rel. error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 320 | 2.0235226 | 2.0943951 | 0.657275 | 3.4×10⁻² |
| 3 | 1280 | 2.0763703 | 2.0943951 | 0.664300 | 8.6×10⁻³ |
| 4 | 5120 | 2.0898693 | 2.0943951 | 0.666077 | 2.2×10⁻³ |
| 5 | 20480 | 2.0932622 | 2.0943951 | 0.666508 | 5.4×10⁻⁴ |

Error falls by ~4× per subdivision, second order in the mesh edge length, as it
should for piecewise-linear interpolation of a smooth field.

This is the strongest single check of the interpolation-and-integration path,
and it is also the canonical application: supply a radiance field and this
integral is a radiant intensity. FlatLand attaches no physical meaning to it —
it is the area integral of a scalar field, and what that scalar means is the
caller's business.

### Convergence with pixel size

Unit cube, generic view direction (0.3, 0.9, −0.31), exact answer 1.51295314.

| h | FlatLand | Abs. error |
| ---: | ---: | ---: |
| 0.02 | 1.52160000 | 8.6×10⁻³ |
| 0.01 | 1.51190000 | 1.1×10⁻³ |
| 0.005 | 1.51192500 | 1.0×10⁻³ |
| 0.0025 | 1.51183750 | 1.1×10⁻³ |
| 0.00125 | 1.51305000 | 9.7×10⁻⁵ |
| 0.000625 | 1.51304844 | 9.5×10⁻⁵ |

Fitted order **h^1.22**, consistent with first-order convergence from a boundary
band one pixel wide.

The error is deliberately shown as **non-monotonic**, because it genuinely is: a
silhouette edge can land favourably on the pixel grid at one resolution and
badly at the next, so halving `h` does not always halve the error. Anyone
choosing a resolution should expect that. The practical guidance is unchanged —
sweep `-r` and watch the answer settle, rather than trusting a single value.

### Invariances

Properties that must hold *identically*, not approximately. These are the ones
whose violation produced silently wrong answers in earlier versions.

| Property | Result |
| :--- | :--- |
| Translating the mesh by 5×10⁶ leaves the area unchanged (float) | exact, 0.00×10⁰ |
| Translating the mesh by 5×10⁶ leaves the area unchanged (double) | exact, 0.00×10⁰ |
| Culled and `--no-cull` agree on area for a closed mesh | exact, 0.00×10⁰ |
| Culled and `--no-cull` agree on the field mean | exact, 0.00×10⁰ |
| Scaling by 3 scales the area by 9 | 1.3×10⁻¹⁶ (machine precision) |

Before the fixes in this release, the translation case reported an area **24.9%
high** in float, and the culling cases disagreed on the field mean because the
depth test resolved to the far surface.

---

## Tolerances

The suite's tolerances are set from these measurements with roughly an order of
magnitude of headroom, tightest where the physics is exact:

| Study | Tolerance | Rationale |
| :--- | :--- | :--- |
| Exact polyhedra | 1.5–4×10⁻³ | boundary band, scaled to the silhouette perimeter |
| Cauchy | 4–10×10⁻³ | averaging cancels most boundary error |
| Sphere vs. its own mesh | 4×10⁻³ | pure rasterization |
| Lambertian, subdiv 2 | 5×10⁻² | coarse mesh dominates |
| Lambertian, subdiv ≥4 | 1×10⁻² | mesh error has largely gone |
| Invariances | 10⁻⁹ or exact | these are identities, not approximations |

---

## What this does *not* establish

Stated plainly, because a validation page that only lists successes is not much
use:

- **Non-convex geometry has no closed form here.** Every area study uses a
  convex body, because that is where $\tfrac14\sum\lvert n\cdot c\rvert$ applies.
  Self-occlusion is exercised by the cube's field cases (the near face hides the
  far one) and by the sphere, but a mesh with deep concavities is only covered
  by the behavioural suite, not by an analytic result.
- **Open and non-manifold surfaces** are checked only through the flat-sheet
  field cases with `--no-cull`.
- **Float precision is under-represented.** Most studies run in double to isolate
  the rasterizer's own error from floating-point noise. The float path is
  covered separately by the invariance and conditioning tests in
  `tests/test_geometry.sh`.
- **These are accuracy studies, not performance ones.** Nothing here measures
  throughput.
- **The mesh generators are shared** between the expectations and the inputs, so
  a bug in `icosphere()` itself would affect both. That is mitigated by
  `self_check.py` pinning the generated shapes against independently-known
  constants (surface areas, edge lengths, convergence limits), but it is not the
  same as an independently sourced mesh.

## Files

| Path | Purpose |
| :--- | :--- |
| `validation/cases.py` | mesh generators and the closed-form expectations |
| `validation/self_check.py` | checks the formulas themselves, before any FlatLand run |
| `validation/run_validation.py` | drives FlatLand and reports the errors |
| `tests/test_validation.sh` | the quick subset, wired into `make test` |

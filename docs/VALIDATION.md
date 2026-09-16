# Validation

FlatLand is checked against closed-form results, not against its own previous
output. Every expected value on this page is derived analytically; nothing is a
stored baseline, so a wrong answer cannot be blessed into a regression fixture.

Reproduce everything here with:

```shell
make validate            # the full study, ~8 s
make test                # the quick subset, as part of the test suite
```

The numbers below were produced on Linux/x86-64 with GCC 13 in double precision.
Field arithmetic will shift in the last digit or two on other toolchains; the
tolerances are set well clear of that.

**Coverage, though, is not allowed to drift at all.** Which pixels a triangle
claims is a discrete decision, so a one-ULP difference there is not a small error
but a whole pixel. The rasterizer decides coverage from three edge functions and
relies on `edge(a,b,p) == -edge(b,a,p)` holding exactly, so that a pixel centre
lying on the shared edge of two triangles is claimed by both rather than by
neither. Fusing a multiply and an add into an FMA breaks that antisymmetry and
opens one-pixel cracks along shared edges — so FlatLand builds with
`-ffp-contract=off` on every compiler that accepts it, and
`tests/test_geometry.sh` pins exact pixel counts across a resolution sweep to
keep it that way. Without the flag, a subdivided unit cube seen head-on covers
98 pixels instead of 100 at `-r 0.1`, on arm64 or on x86-64 built with `-mfma`.

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


### Reading the integral as a radiant intensity

Three of the studies below are phrased as blackbody radiation. That is not a
feature of the tool — it is one reading of the same generic integral, and it is
worth spelling out because it is the reading most users arrive with.

FlatLand integrates over the *projected* area, and projected area is
$dA_\perp = \cos\theta\,dA$ with $\theta$ the angle between the surface normal
and the viewer. So for a field $f = L$, a radiance,

$$I(\hat n) = \int L\,dA_\perp = \int_{\text{visible}} L\cos\theta\,dA
\qquad [\mathrm{W/sr}]$$

which is the definition of radiant intensity. A blackbody is a Lambertian
emitter, so its radiance is isotropic and fixed by temperature alone:

$$M = \sigma T^4 \quad [\mathrm{W/m^2}], \qquad L = M/\pi \quad [\mathrm{W/m^2\,sr}]$$

The $\pi$ is the projected solid angle of a hemisphere,
$\int\cos\theta\,d\Omega = \pi$. Two consequences make good test cases: the
intensity of an *isothermal* body is its projected area times $\sigma T^4/\pi$
from every direction, and integrating that over all directions must return the
Stefan–Boltzmann law.

$\sigma$ itself is not a fitted number. The 2019 SI fixes $h$, $k$ and $c$
exactly, so $\sigma = 2\pi^5k^4/15h^3c^2 =
5.670374419\!\times\!10^{-8}\ \mathrm{W\,m^{-2}K^{-4}}$ is exact, and
`self_check.py` rebuilds it from those three constants rather than trusting the
literal in `cases.py`.

The formulas themselves are checked before anything else runs
(`validation/self_check.py`, 67 identities) against values known from elsewhere —
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

### Blackbody radiation — closed forms with a physical reading

Four cases, ordered by how much they ask of the tool. All run at 1200 K in
double precision. See [Reading the integral as a radiant
intensity](#reading-the-integral-as-a-radiant-intensity) for the radiometry.

**1. Isothermal sphere.** The radiance is constant, so $A_\perp = \pi r^2$ gives

$$I = \frac{\sigma T^4}{\pi}\,\pi r^2 = \sigma T^4 r^2
\qquad\text{from every direction}$$

Twelve Fibonacci directions on an icosphere(5), exact value 117580.884 W/sr:

| Quantity | Spread over 12 directions | Rel. error |
| :--- | ---: | ---: |
| $I$ [W/sr] | 117543.84 – 117548.00 | 2.8–3.2×10⁻⁴ |
| $L$ [W/m²/sr] | 37427.15779 (identical to 11 digits) | 1.4×10⁻¹¹ |

The two rows differ by seven orders of magnitude, and the reason is the whole
argument of this page. $I$ inherits the icosphere's projected-area deficit — the
mesh is inscribed, so it genuinely projects to slightly less than $\pi r^2$.
$L$ is the *mean*, and a constant field must interpolate to exactly that constant
at every covered pixel regardless of what shape those pixels cover. So the mean
measures interpolation alone, and it is exact.

**2. Stefan–Boltzmann, recovered from projected areas.** Integrating the
intensity over all directions and applying Cauchy's identity:

$$\oint I\,d\Omega = \frac{\sigma T^4}{\pi}\oint A_\perp\,d\Omega
= \frac{\sigma T^4}{\pi}\,(4\pi)\frac{S}{4} = \sigma T^4 S$$

The left-hand side never uses the surface area; getting $\sigma T^4S$ back out is
the check. 800 directions, resolution 2×10⁻³:

| Shape | Triangles | $4\pi\langle I\rangle$ [W] | $\sigma T^4 S$ [W] | Rel. error |
| :--- | ---: | ---: | ---: | ---: |
| cube | 12 | 705497.75 | 705485.30 | 1.8×10⁻⁵ |
| octahedron | 8 | 814633.47 | 814624.26 | 1.1×10⁻⁵ |
| icosphere(3) | 1280 | 1470538.88 | 1470524.47 | 9.8×10⁻⁶ |

The residual is direction *sampling*, not FlatLand — and it is measurably so.
Halving the pixel size to 10⁻³ moves the worst case from 1.8×10⁻⁵ to 1.1×10⁻⁵
and costs four times as much; the number of directions is what sets it, which is
the same quantity the Cauchy study measures. This case therefore runs at a
coarser resolution than the rest of the study on purpose.

**3. Graded $T^4$ — affine in position, no terminator.** Take a sphere with

$$T(p)^4 = T_{\max}^4\,\frac{1 + \hat p\cdot\hat s}{2}$$

hot at the pole facing $\hat s$, falling to absolute zero at the antipode,
non-negative everywhere. The radiance is then *linear in the vertex coordinates*,
so barycentric interpolation reproduces it exactly and no field-representation
error enters at all. Two standard hemisphere integrals,
$\int\cos\theta\,d\Omega = \pi$ and
$\int(\hat p\cdot\hat s)\cos\theta\,d\Omega = \tfrac{2\pi}{3}\cos\alpha$, give

$$I(\alpha) = \frac{\sigma T_{\max}^4 r^2}{2}\left[1 + \tfrac{2}{3}\cos\alpha\right]$$

with $\alpha$ the phase angle between the source and the observer. icosphere(5),
20480 triangles, resolution 10⁻³:

| Phase | $I$ FlatLand | Exact for the *mesh* | Exact for the *sphere* | Raster err. | vs. sphere |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0° | 97945.042 | 97944.156 | 97984.070 | 9.1×10⁻⁶ | 4.0×10⁻⁴ |
| 30° | 92697.474 | 92696.801 | 92733.119 | 7.3×10⁻⁶ | 3.8×10⁻⁴ |
| 60° | 78359.804 | 78359.429 | 78387.256 | 4.8×10⁻⁶ | 3.5×10⁻⁴ |
| 90° | 58772.614 | 58771.724 | 58790.442 | 1.5×10⁻⁵ | 3.0×10⁻⁴ |
| 120° | 39187.382 | 39186.996 | 39193.628 | 9.9×10⁻⁶ | 1.6×10⁻⁴ |
| 150° | 24848.832 | 24848.157 | 24847.764 | 2.7×10⁻⁵ | 4.3×10⁻⁵ |

The middle column is the exact integral over the polyhedron FlatLand was
actually handed (`cases.field_integral_convex`), so the "raster err." column is
rasterization alone — ~10⁻⁵ — while "vs. sphere" adds the mesh's own fidelity.
Same separation as the sphere area study, same conclusion: the tool reproduces
its input mesh to five digits, and the rest belongs to the mesh.

**4. Radiative equilibrium — a real terminator.** A sphere in instantaneous
equilibrium with a distant source balances absorbed flux $\propto\cos\psi$
against $\sigma T^4$, giving the classic subsolar law

$$T(\psi) = T_{\text{sub}}\cos^{1/4}\psi \quad\text{on the lit side},\qquad 0\ \text{beyond it}$$

so $\sigma T^4 = \sigma T_{\text{sub}}^4\max(\cos\psi, 0)$. The disc-integrated
result is the Lambert-sphere phase function — Russell's 1916 planetary
photometry result:

$$I(\alpha) = \tfrac{2}{3}\sigma T_{\text{sub}}^4 r^2\,\Phi(\alpha),
\qquad \Phi(\alpha) = \frac{\sin\alpha + (\pi-\alpha)\cos\alpha}{\pi}$$

with $\Phi(0)=1$, $\Phi(\pi/2)=1/\pi$, $\Phi(\pi)=0$. Unlike case 3 this has a
**kink** at the terminator, which per-vertex linear interpolation cannot
represent, so the mesh column should be worse. Same mesh and resolution:

| Phase | $I$ FlatLand | Exact for the *mesh* | Exact for the *sphere* | Raster err. | vs. sphere |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0° | 78344.855 | 78344.865 | 78387.256 | 1.2×10⁻⁷ | 5.4×10⁻⁴ |
| 30° | 69012.250 | 69011.912 | 69046.848 | 4.9×10⁻⁶ | 5.0×10⁻⁴ |
| 60° | 47716.559 | 47716.273 | 47737.665 | 6.0×10⁻⁶ | 4.4×10⁻⁴ |
| 90° | 24940.786 | 24940.442 | 24951.439 | 1.4×10⁻⁵ | 4.3×10⁻⁴ |
| 120° | 8544.137 | 8543.840 | 8544.037 | 3.5×10⁻⁵ | 1.2×10⁻⁵ |
| 150° | 1163.591 | 1163.269 | 1161.493 | 2.8×10⁻⁴ | 1.8×10⁻³ |

It is, and in the expected place: at 150° only a thin lit crescent is visible,
the terminator runs through most of it, and the error against the true sphere
grows to 1.8×10⁻³ while the rasterization column stays at 10⁻⁴. The 1.2×10⁻⁵ at
120° is a **coincidence** — the mesh's area deficit and the terminator error have
opposite signs and cross there. That is exactly why the check is a sweep and not
a single angle; agreement at one point is not evidence.

Against the smooth sphere both cases converge fourfold per subdivision level at
low phase — second order in the mesh edge length, as piecewise-linear
interpolation requires. Measured at 0° phase for subdivisions 2–5:

| Subdiv | Triangles | Graded, vs. sphere | Equilibrium, vs. sphere |
| ---: | ---: | ---: | ---: |
| 2 | 320 | 2.56×10⁻² | 3.38×10⁻² |
| 3 | 1280 | 6.49×10⁻³ (×3.9) | 8.61×10⁻³ (×3.9) |
| 4 | 5120 | 1.63×10⁻³ (×4.0) | 2.16×10⁻³ (×4.0) |
| 5 | 20480 | 4.07×10⁻⁴ (×4.0) | 5.41×10⁻⁴ (×4.0) |

The graded case holds that rate at every phase angle. The equilibrium case does
not: at 120° it refines by only ×2.3 between subdivisions 3 and 4, because the
terminator error is a different and slower-converging term that the smooth case
does not have. That is the cost of the kink, and it is visible in the data rather
than argued for.

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
| Blackbody, mean radiance | 10⁻⁹ | a constant field is interpolated exactly |
| Blackbody vs. the mesh's own integral | 1.5×10⁻³ | pure rasterization |
| Blackbody vs. the smooth sphere | 1.5–6×10⁻³ | mesh fidelity dominates |
| Stefan–Boltzmann closure | 10⁻⁴ | direction sampling, scaled to the count |
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
- **The blackbody cases cover emission only, and only spectrally integrated.**
  $\sigma T^4$ is the *total* exitance over all wavelengths. A band-limited
  intensity needs the Planck fraction $F(\lambda T)$, which has no elementary
  closed form — it is a rapidly convergent series, not an identity — so no study
  here pins one. Nothing here covers absorption, reflection, or any second
  bounce either: FlatLand resolves one visible surface per pixel and integrates
  over it, which is what these cases test and all they test.
- **The quick tier cannot catch a wrong formula.** Its mesh tolerances are sized
  for a coarse icosphere (1–3×10⁻²), so a percent-level error in a closed form
  would pass it. The full study does catch one, and `self_check.py` — which runs
  first, in both tiers — catches it decisively. That layering is the defence, not
  the study tolerances.
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

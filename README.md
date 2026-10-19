# Core entropy of the real quadratic family

An interactive plot of the core entropy `h(λ)` of

```
g_λ(x) = λ x (1 − x),   2 ≤ λ ≤ 4,   x ∈ [0,1]
```

computed by a small dependency-free Julia module and drawn by two front ends:
an interactive HTML page and a Makie window.

Files:

| file | what |
|---|---|
| `core_entropy.jl` | front end: runs the sweep, writes the CSV and the HTML |
| `src/CoreEntropy.jl` | the computational core — kneading theory, entropy solver, grid, CSV, shared helpers |
| `src/checks.jl` | the validation suite (`--check`) |
| `src/HTMLPage.jl` | the Plotly page |
| `core_entropy_makie.jl` | front end: the Makie window |
| `core_entropy.html` | the interactive plot (plotly.js embedded, works offline, ~5 MB) |
| `core_entropy_data.csv` | `lambda,c,h` — 28 696 rows, 12 decimals |
| `preview.png` | static preview |
| `preview_bifurcation.png` | static preview with the bifurcation overlay |
| `core_entropy_notes.pdf` | ten-page write-up of the mathematics and the algorithm, with figures (source: `core_entropy_notes.tex`) |
| `figs/` | the six figures of the write-up (`make_figures.py` rebuilds them from the CSV) |

Both front ends compute with the same `src/CoreEntropy.jl`; neither of them contains any
mathematics of its own.

## The function

`g_λ` is affinely conjugate to `z ↦ z² + c` with

```
c = λ/2 − λ²/4  ∈ [−2, 0]        (λ = 1 + sqrt(1 − 4c))
```

For real `c` the Hubbard tree of the quadratic is the real interval spanned by the
postcritical set, and the **core entropy** is the topological entropy of the map
restricted to it; for `λ ≥ 2` this equals the topological entropy of `g_λ` on `[0,1]`.
Tiozzo (*Topological entropy of quadratic polynomials and dimension of sections of the
Mandelbrot set*, Adv. Math. 273, 2015) proved that `h` is continuous in `c`; here it is
computed rather than quoted.

`h` is a devil's staircase: nondecreasing, constant on every hyperbolic component, and
it is zero on the whole period-doubling cascade, i.e. for `2 ≤ λ ≤ λ∞ ≈ 3.5699457`
(`c ∈ [−1.4011552, 0]`). Landmarks:

| λ | c | h |
|---|---|---|
| 2 | 0 | 0 (critical point fixed) |
| 3 | −3/4 | 0 (parabolic fixed point) |
| 1+√5 ≈ 3.2360679 | −1 | 0 (basilica centre) |
| λ∞ ≈ 3.5699457 | ≈ −1.4011552 | 0 (Feigenbaum accumulation) |
| ≈ 3.8318741 | ≈ −1.7548777 | log φ ≈ 0.4812118 (period-3 centre) |
| 4 | −2 | log 2 ≈ 0.6931472 (Chebyshev map) |

The step at `h = log φ` is unusually wide. The period-3 hyperbolic component itself is
`[3.8284, 3.8415]` (saddle-node at `1+√8`, period doubling at `3.8414990`), but `h` stays
exactly `log φ` out to `λ ≈ 3.8568`, i.e. over `c ∈ [−1.789, −1.75]` — the real extent of
the period-3 copy of the Mandelbrot set — and only then rises again. At `λ = 3.845` the
kneading determinant is exactly `(1 − t − t²)·Q(t)/(1 − t¹²)` with `Q` of degree 9, so the
entropy is exactly `log φ`; a 60-digit evaluation confirms `D(1/φ) = 0` at `λ = 3.85,
3.855, 3.8565` as well.

## How h is computed

Let `x_j = g_λ^j(1/2)` be the critical orbit and

```
s_j = sign g_λ'(x_j) = +1 if x_j < 1/2,  −1 if x_j > 1/2
```

(a return `x_j = 1/2` is read as the `0`-symbol, `s_j = −1`). With
`ε_k = s_1 s_2 ⋯ s_k`, the Milnor–Thurston kneading determinant is

```
D(t) = Σ_{k≥0} ε_k t^k ,
```

and the entropy is `h = log(1/t₀)`, where `t₀` is the smallest zero of `D` in `(0,1)`:
`1/t₀` is the growth rate of the number of laps of `g^n`. If `D` has no zero in `(0,1)`
then `h = 0`. The truncation error is bounded by the tail, since `|ε_k| = 1`:

```
|D(t) − D_N(t)| ≤ t^(N+1)/(1−t) .
```

`core_entropy` therefore accepts a zero only once `(1−t₀)N ≥ 25` (so the tail is below
`e^{-25}`), doubling `N` from 800 up to 51 200 when `h` is small; `h` is reported as 0,
resolved down to `h ≈ 0.008`, when no zero exists at that resolution. Over a grid the
zero from the previous parameter is used as the initial bracket, which is what makes the
28 696-parameter sweep take about 1.5 s.

## Running it

```sh
julia core_entropy.jl          # recompute + rewrite core_entropy_data.csv and core_entropy.html
julia core_entropy.jl --check  # validation, no files written
julia core_entropy.jl --cdn    # smaller HTML that loads plotly.js from the CDN
```

Only Julia's standard library is used (`Printf`, `Downloads`), tested with Julia 1.11.3.
The first run downloads plotly.js 2.35.2 into `~/.cache/core_entropy/` and embeds it in
the HTML, so the result is a single self-contained file; later runs reuse the cache. If
the download fails the HTML falls back to a CDN `<script>` tag.

## The interactive page

* the curve is plotted in units of `log 2`, so that its maximum `h(4) = log 2` sits at `1` and
  the axis reads as a fraction of the maximal entropy; hover for `λ`, `c`, `h/log 2` and the
  raw `h`; scroll/drag to zoom, the bottom range slider zooms too, and the top axis carries
  the `c` scale;
* click anywhere on the curve, drag the slider, or use the landmark links to choose `λ`;
* two checkboxes control the view.  **entropy curve** (on by default) hides the curve,
  its marker and the left axis, leaving the `λ` axis and whatever else is shown;
  **bifurcation diagram** adds the attracting set of `g_λ` as a greyscale density — 1000
  parameter columns × 400 bins in `x`, 1500 transient and 1200 plotted iterates each,
  log-scaled — read on the right-hand axis (`x ∈ [0,1]`).  With the curve off and the
  diagram on, the page is an ordinary bifurcation plot over `[2,4]`; with both on, the
  plateaus of `h` sit visibly on the corresponding windows.  The diagram is computed in
  the browser for the visible `λ`-range, recomputed on zoom, and costs nothing when off;
* the lower panel shows `g_λ` on `[0,1]`, the critical orbit, and the core (shaded), and
  reports the itinerary of the critical orbit (`L`/`R` of the critical point) and its
  eventual period.

## Validation

`julia core_entropy.jl --check` reports:

```
h(2)          = 0                        exact, critical point fixed
h(3)          = 0                        exact, parabolic fixed point
h(1+√5)       = 0                        exact, basilica centre
h(3.8318741)  = 0.481211825060 = log φ   |err| < 1e-16
h(4)          = 0.693147180560 = log 2   |err| < 1e-16
```

plus an independent check that the number of laps of `g^n` (preimages of the critical
point, computed from the dynamics alone) grows like `e^{hn}`, and monotonicity of `h` on
a 4001-point sweep. The `log φ` plateau was checked separately: at `λ = 3.845` the sign
sequence is exactly 12-periodic and the numerator of `D` is exactly divisible by
`t² + t − 1`, so the zero is exactly `1/φ`; at `λ = 3.85, 3.855, 3.8565` a 60-digit
evaluation gives `D(1/φ) = 0` to the precision of `φ` itself.

## A native window

`core_entropy_makie.jl` is the Makie counterpart of the HTML page: the entropy curve
with a `λ` slider that drives a panel showing `g_λ` on `[0,1]` and its critical orbit.
It includes only `src/CoreEntropy.jl` — not the HTML machinery — and contains nothing but
the drawing, so the two front ends cannot drift apart. It needs no CSV.

```sh
julia -e 'using Pkg; Pkg.add(["Makie", "GLMakie"])'   # once
julia core_entropy_makie.jl
```

It has **not been run**: the machine this was written on has neither a display nor
Makie, so the plotting calls are untested (the entropy code it calls is the tested one;
the helpers were checked against the real grid). It assumes a recent Makie (`Slider`
with `startvalue`, `Axis` with `limits`). `core_entropy_data.csv` has the columns
`lambda,c,h` if you would rather build your own figure.

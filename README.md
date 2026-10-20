# Core entropy of the real quadratic family

An interactive plot of the core entropy `h(λ)` of

```
g_λ(x) = λ x (1 − x),   2 ≤ λ ≤ 4,   x ∈ [0,1]
```

computed by a small dependency-free Julia module and drawn by two front ends:
an interactive HTML page and a Makie window.

Files.  The generated ones are marked *built*: they are not in the repository, and `make`
produces them from the sources.

| file | what |
|---|---|
| `src/CoreEntropy.jl` | the computational core — kneading theory, entropy solver, grid, CSV, shared helpers |
| `src/checks.jl` | the validation suite (`--check`) |
| `src/HTMLPage.jl` | the Plotly page |
| `core_entropy.jl` | front end: runs the sweep, writes the CSV and the HTML |
| `core_entropy.py` | the same sweep in Python (numpy); agrees with the Julia one to 5e-11 |
| `tent_map.py` | Python front end: entropy, bifurcation diagram and tent map on one page |
| `tent_map_app.py` | Python front end: the same, in a native matplotlib window |
| `core_entropy_makie.jl` | front end: the Makie window |
| `core_entropy_notes.tex` | the write-up |
| `bands.tex` | the companion note: the bands and their renormalization |
| `figs/make_figures.py`, `figs/make_previews.py` | the figure and preview sources |
| `Makefile` | builds everything below |
| `core_entropy_data.csv` | *built* — `lambda,c,h`, 28 696 rows, 12 decimals |
| `core_entropy.html` | *built* — the interactive entropy page (plotly.js embedded, ~5 MB) |
| `tent_map.html` | *built* — the integrated page (~5 MB) |
| `figs/fig_*.pdf` | *built* — the seven figures of the write-up and the three of the companion note |
| `core_entropy_notes.pdf` | *built* — the fourteen-page write-up |
| `bands.pdf` | *built* — the five-page companion note |
| `preview*.png` | *built* — static previews |

## Building

```sh
make            # everything (about 30 s)
make note       # just the write-up, from the figures and the .tex
make bands      # just the companion note on bands and renormalization
make app        # the native matplotlib window (needs a display)
make check      # validation of the Julia and the Python implementations
make clean      # remove every generated file
```

`make` needs Julia (`core_entropy.jl`), Python 3 with numpy and matplotlib, `latexmk` with
the usual TeX packages, and `pdftoppm` from poppler for one preview; `JULIA=`, `PYTHON=` and
`LATEXMK=` override the tools.

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
Its continuity in the parameter is Milnor and Thurston's (Lemma 12.3: the growth number is
continuous in the map at fixed lap number). For real parameters Tiozzo proved the formula
relating the entropy to the dimensions of sections of the Mandelbrot set (*Topological entropy
of quadratic polynomials and dimension of sections of the Mandelbrot set*, Adv. Math. 273,
2015), and he later proved continuity of the core entropy for all quadratic polynomials, as a
function of the external angle (*Continuity of core entropy of quadratic polynomials*, Invent.
Math. 203, 2016). Here `h` is computed rather than quoted.

`h` is a devil's staircase: nondecreasing, constant on every hyperbolic component, and
it is zero on the whole period-doubling cascade, i.e. for `2 ≤ λ ≤ λ∞ ≈ 3.5699457`
(`c ∈ [−1.4011552, 0]`). Landmarks:

| λ | c | h |
|---|---|---|
| 2 | 0 | 0 (critical point fixed) |
| 3 | −3/4 | 0 (parabolic fixed point) |
| 1+√5 ≈ 3.2360679 | −1 | 0 (basilica centre) |
| λ∞ ≈ 3.5699457 | ≈ −1.4011552 | 0 (Feigenbaum accumulation) |
| ≈ 3.8318741 | ≈ −1.7548777 | log τ ≈ 0.4812118 (period-3 centre) |
| 4 | −2 | log 2 ≈ 0.6931472 (Chebyshev map) |

The step at `h = log τ` is unusually wide. The period-3 hyperbolic component itself is
`[3.8284, 3.8415]` (saddle-node at `1+√8`, period doubling at `3.8414990`), but `h` stays
exactly `log τ` out to `λ ≈ 3.8568`, i.e. over `c ∈ [−1.789, −1.75]` — the real extent of
the period-3 copy of the Mandelbrot set — and only then rises again. At `λ = 3.845` the
kneading determinant is exactly `(1 − t − t²)·Q(t)/(1 − t¹²)` with `Q` of degree 9, so the
entropy is exactly `log τ`; a 60-digit evaluation confirms `D(1/τ) = 0` at `λ = 3.85,
3.855, 3.8565` as well.  The same flat step appears at the period-5 window,
`λ ≈ 3.9056–3.9068`, where `h = log s₅` with `s₅⁴ − s₅³ − s₅² − s₅ + 1 = 0`
(`s₅ = 1.7220838…`), and at the period-4 window, `λ ≈ 3.9601–3.9616`, where
`h = log s₄` with `s₄³ − s₄² − s₄ − 1 = 0` (`s₄ = 1.8392868…`, the tribonacci constant).

## How h is computed

Let `x_j = g_λ^j(1/2)` be the critical orbit and

```
s_j = sign g_λ'(x_j) = +1 if x_j < 1/2,  −1 if x_j > 1/2
```

(a return `x_j = 1/2` is read as the `1`-symbol, `s_j = −1`). With
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

`core_entropy` therefore accepts a zero only once `hN ≥ 25` (so the tail is below
`e^{-25}/h`), raising `N` adaptively from 800 up to 51 200 when `h` is small; `h` is reported as 0,
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

* the curve is plotted in nats, with the axis running from 0 to just above `log 2`, so that
  the maximum `h(4) = log 2` sits at the top of the frame; hover for `λ`, `c`, `h` and
  `h/log 2`; scroll/drag to zoom, the bottom range slider zooms too, and the top axis carries
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
h(3.8318741)  = 0.481211825060 = log τ   |err| < 1e-16
h(4)          = 0.693147180560 = log 2   |err| < 1e-16
```

plus an independent check that the number of laps of `g^n` (preimages of the critical
point, computed from the dynamics alone) grows like `e^{hn}`, and monotonicity of `h` on
a 4001-point sweep. The `log τ` plateau was checked separately: at `λ = 3.845` the sign
sequence is exactly 12-periodic and the numerator of `D` is exactly divisible by
`t² + t − 1`, so the zero is exactly `1/τ`; at `λ = 3.85, 3.855, 3.8565` a 60-digit
evaluation gives `D(1/τ) = 0` to the precision of `τ` itself.

## The integrated page

`tent_map.py` (Python, numpy only) writes `tent_map.html`, which puts the whole story on one
`λ` axis:

- the core entropy `h(λ)`, computed by `core_entropy.py`, with the **bifurcation diagram**
  of the same parameters behind it as a greyscale density on the right-hand axis — togglable,
  and recomputed for the visible range when you zoom;
- the cobweb of the critical orbit for `g_λ` and for the tent map `f_s` of slope
  `s = e^{h(λ)} = 1/t₀`, with their itineraries — or, with the *plot the iterates* toggle,
  the graphs of `g_λ^n` and `f_s^n` themselves for `n = 1…10`.  The `n`-th iterate has
  `ℓ_n = 1 +` its turning points monotone pieces, so watching the count climb — 4, 14, 38,
  90 at `λ = 3.7`, the same for both maps — is the most concrete picture of the quantity
  whose growth rate is the entropy;
- the map `A` of Theorem 7.4 of Milnor–Thurston, which satisfies `f_s(A(x)) = A(g_λ(x))`:
  `g_λ` is topologically semi-conjugate to the piecewise-linear model `f_s`.  The arrow does
  not reverse — a tent map has no attracting cycle, while `g_λ` may have one — and `A` need
  not be injective: where it collapses an interval, the two itineraries can differ.  The
  write-up's §4.3 gives the exact criterion for `A` to be a conjugacy.

Everything follows the slider.  `core_entropy.py` ports `src/CoreEntropy.jl` to numpy
(a vectorised sweep over the parameters, with the scalar adaptive-order fallback), and its
values agree with the Julia ones to `5e-11` on all 28 696 parameters; `tent_map.py` reuses
`core_entropy_data.csv` when it matches the grid and recomputes it otherwise.

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

`tent_map_app.py` is the same window in Python and matplotlib, so it needs neither Julia nor
a browser: the entropy curve with a `λ` slider, `g_λ` on `[0,1]` with its critical orbit, and
the tent map `f_s` of slope `s = e^{h(λ)}` drawn in the same panel.

```sh
python3 tent_map_app.py                              # open the window
python3 tent_map_app.py --exact                      # h computed on the fly
python3 tent_map_app.py --save frame.png --lam 3.7   # one frame, no display
```

`h` exists only at the grid parameters, so by default the slider snaps to the nearest one and
displays which it used — near the merging parameter the grid spacing alone moves `h` by
`3·10⁻⁴`, and a readout at an off-grid `λ` would misreport it. `--exact` computes `h` at the
slider's own `λ` instead and prints how long that took: about 10 ms above `λ = 3.6` and about
36 ms below, where the solver must quadruple `N` before it can report zero. A full redraw
costs 25 ms on top, so a drag runs at 16–20 updates a second in either mode. Like the Makie window this has **not been driven
by a real GUI backend** here; `--save` exercises the drawing, and the callbacks were run
headlessly across `λ` from 2 to 4.

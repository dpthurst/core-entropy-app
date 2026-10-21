#!/usr/bin/env python3
"""
tent_map_app.py -- core entropy, bifurcation diagram and the semi-conjugate tent
map of the logistic family, in a native window.

The matplotlib counterpart of tent_map.html (and of core_entropy_makie.jl on the
Julia side): drag the slider to move lambda, and zoom or pan the lambda axis to
look at a range.  No browser, no plotly, no Julia -- the only requirement beyond
numpy is matplotlib.

    python3 tent_map_app.py                              # open the window
    python3 tent_map_app.py --lam 3.7                    # start elsewhere
    python3 tent_map_app.py --exact                      # h computed on the fly
    python3 tent_map_app.py --save frame.png --lam 3.7   # one frame, headless
    python3 tent_map_app.py --save zoom.png --xlim 3.56 3.60

Left: h(lambda), with the bifurcation diagram of the same parameters behind it as
a greyscale density read on the right-hand axis.  Zooming or panning that axis
recomputes the diagram for the visible lambda range, which is what makes the small
scales readable, and rescales the slider beneath it to the same window -- so after
zooming, the slider's two ends are the ends of the visible range and every pixel of
it is spent inside it.  Zoom and pan act on lambda alone: the toolbar's rectangle
applies only its width -- the band is drawn at the full height to say so -- and a pan
ignores the vertical component, so the entropy axis stays at [-0.02, 0.70] for
comparison.  The window opens with the zoom tool selected.  Every visited cell of the diagram is drawn with at least half ink, so a
window containing an attracting cycle -- where the counts run from 1 to 700 -- reads as
well as the full range.  Right: g_lambda and the tent map f_s of slope s = exp(h(lambda)),
each with the cobweb of the critical orbit.

Right, three panels above the slider: g_lambda and the tent map f_s, each with the
cobweb of the critical orbit -- or, with the "iterates" toggle and its n, the graphs of
the n-th iterates themselves, titled with their lap counts, which the two maps share
(the + and - keys, or the up and down arrows, step n one at a time);
and Milnor-Thurston's semi-conjugacy A, with the images of the first seven critical
points on its graph.  A sits directly above g and the tent map sits to its right with its
argument on the vertical axis, so the identity f_s(A(x)) = A(g(x)) can be read off the
picture: a vertical segment from (x, g(x)) up to (x, A(x)), then a horizontal one at that
height to the tent graph, which meets it at f_s(A(x)) -- the value A takes at the image
point.  In the cobweb view there is one such connection at each of the first few points of
the orbit, and the length of the cobweb is set by its own slider; showing one per point
would be a mess, and a few is enough to see A carrying the g-orbit onto the f_s-orbit.  In
the iterates view there is a single one, at a generic x, and it is drawn only at n = 1,
where the panels really do show g and f_s; nothing is drawn with A switched off.  Two
checkboxes switch the iterates view (on at start, at n = 1) and the A panel (off at
start); with A off a slider event costs 48 ms rather than 60.  A red vertical line marks
the current lambda.

By default h is read off the 28 696-parameter grid, which is instant; with --exact
it is computed at the slider's own lambda by core_entropy.py, at about 10 ms above
lambda = 3.6 and 36 ms below, where the solver must quadruple N before it can
report zero.  Results are cached, so revisiting a value is free.  A recomputed
bifurcation diagram costs about 10-15 ms, and is built once per zoom rather than
per mouse movement.

A slider event costs one full redraw, about 60 ms here -- the four panels cost 7-13 ms
each, the A panel 9 ms of it, and that is the floor for a pyplot window this size
without blitting; switching A off brings it back to 48 ms.  The
Slider's own drawing is switched off (sl.drawon = False) because set_val draws
*before* it notifies its observers, so leaving it on renders the figure twice, and
the first time with the panels still stale.

NOT TESTED INTERACTIVELY: written in an environment with no display, so the window
and the mouse have never been driven by a real GUI backend.  The drawing code is
exercised by --save, and the slider and zoom callbacks were driven headlessly.
"""
import argparse
import os
import time

import numpy as np

import core_entropy as ce

HERE = os.path.dirname(os.path.abspath(__file__))
LOG2 = np.log(2.0)
BLUE, GREEN, RED, ORANGE = "#1c5f9e", "#2f6b4f", "#c1272d", "#e08a1e"
LAM_LO, LAM_HI = 2.0, 4.0

# Tuning, kept here rather than in the call signatures: the diagram's grid is about
# twice the panel's pixels (more is invisible, and the image is redrawn on every slider
# event), and A is truncated at 64 foldings, which leaves a residual of 1e-11 or better.
BIF_NLAM, BIF_NROWS, BIF_TRANSIENT, BIF_KEEP = 800, 420, 1200, 700
BIF_GAMMA, BIF_FLOOR = 0.35, 0.6  # ink = floor + (1-floor)*(count/top)**gamma
A_STEPS, A_POINTS = 64, 241
COBWEB_STEPS, ORBIT_DOTS = 60, 7
Y_PIN = (-0.02, 0.70)    # the entropy axis does not move: only lambda is zoomed


def sweep():
    """h(lambda) on the grid, reused from core_entropy_data.csv when it matches."""
    path = os.path.join(HERE, "core_entropy_data.csv")
    lams = ce.lambda_grid()
    if os.path.isfile(path):
        d = np.loadtxt(path, delimiter=",", skiprows=1)
        if d.shape[0] == len(lams) and np.allclose(d[:, 0], lams, atol=1e-9):
            print(f"reusing the sweep in {os.path.basename(path)}")
            return d[:, 0], d[:, 2]
    print("computing the sweep in Python (about 40 s) ...")
    hs = ce.compute_curve(lams, verbose=False)
    return lams, hs


def tent(x, s):
    """The tent map of slope s: s*x up to 1/2, then s*(1-x)."""
    return np.where(x <= 0.5, s * x, s * (1.0 - x))


def critical_orbit(lam, n=ORBIT_DOTS):
    """The first n points of the orbit of the critical point under g_lam."""
    out, x = [], 0.5
    for _ in range(n):
        out.append(x)
        x = lam * x * (1.0 - x)
    return out


def cobweb(f, x0, n):
    """The zigzag of the first n iterates: x -> (x,f(x)) -> (f(x),f(x)) -> ..."""
    xs, ys = [], []
    x = x0
    for _ in range(n):
        y = float(f(x))
        xs += [x, x, y]
        ys += [x, y, y]
        x = y
    return xs, ys


def bifurcation(lo, hi, nlam=BIF_NLAM, nrows=BIF_NROWS, transient=BIF_TRANSIENT,
                keep=BIF_KEEP):
    """Occupancy grid of the attracting set over the visible lambda range.

    Vectorised over lambda: every column is iterated at once, which is what makes a
    recomputation on zoom cheap (about 60 ms) rather than a second or two.
    """
    lams = np.linspace(lo, hi, nlam)
    x = np.full(nlam, 0.5)
    for _ in range(transient):
        x = lams * x * (1.0 - x)
    rows = np.empty((keep, nlam), dtype=np.int32)
    for k in range(keep):
        x = lams * x * (1.0 - x)
        rows[k] = np.clip((x * nrows).astype(np.int32), 0, nrows - 1)
    cols = np.broadcast_to(np.arange(nlam, dtype=np.int32), rows.shape)
    z = np.zeros((nrows, nlam))
    np.add.at(z, (rows, cols), 1.0)
    return z


def preimages_g(y, lam):
    """The real x with g_lam(x) = y (one, when y is the critical value)."""
    disc = 1.0 - 4.0 * y / lam
    if disc <= 0.0:
        return []
    r = np.sqrt(disc)
    return [0.5 * (1.0 + r), 0.5 * (1.0 - r)]


def preimages_tent(y, s):
    """The y' with f_s(y') = y.  None when y is above the peak s/2, which is the image."""
    if y > s / 2.0 + 1e-12:
        return []
    return [v for v in (y / s, 1.0 - y / s) if 1e-12 < v < 1.0 - 1e-12]


def lap_count(lam, s, n, which="g"):
    """1 + the number of distinct interior preimages of 1/2 under the first n-1 iterates.

    That is the lap number of the n-th iterate, counted from the tree rather than off a
    grid, so it does not depend on a plotting resolution.
    """
    dedupe = lambda v: list(np.sort(np.unique(np.round(v, 9))))
    S, allp = [0.5], [0.5]
    for _ in range(max(n - 1, 0)):
        nxt = []
        for y in S:
            nxt += preimages_g(y, lam) if which == "g" else preimages_tent(y, s)
        S = dedupe(nxt)
        allp = dedupe(allp + S)
    return len([p for p in allp if 0.0 < p < 1.0]) + 1


def semiconjugacy(x, lam, s, n=A_STEPS):
    """Milnor-Thurston's A: the itinerary of x under g_lam, folded back through f_s.

    A is the limit of these foldings; it satisfies f_s(A(x)) = A(g_lam(x)) up to the
    truncation, and the residual is what the readout reports.  Vectorised over x, so a
    whole graph costs about a millisecond.
    """
    X = np.atleast_1d(np.asarray(x, dtype=float)).copy()
    w = np.empty((n, X.size))
    for k in range(n):
        w[k] = np.where(X < 0.5, 1.0, -1.0)
        X = lam * X * (1.0 - X)
    u = np.full(X.size, 0.5)
    for k in range(n - 1, -1, -1):
        u = np.where(w[k] > 0, u / s, 1.0 - u / s)
    return u


def ink(z, gamma=BIF_GAMMA, floor=BIF_FLOOR):
    """Map an occupancy grid to ink fractions.

    A bifurcation diagram is nearly a set: a cell belongs to the attracting set or it does
    not, and a column at a periodic parameter piles all 700 iterates into three cells while
    its neighbour spreads them over hundreds.  Normalising the counts therefore either
    washes out the sparse columns or flattens the dense ones -- at lambda in [3.82, 3.875]
    the counts run from 1 to 700.  So every visited cell gets at least `floor` ink, and the
    count modulates the rest: the picture is readable at every zoom, and still shows where
    the orbits linger.
    """
    visited = z > 0
    nz = z[visited]
    top = max(1.0, float(np.percentile(nz, 99.5))) if nz.size else 1.0
    out = np.zeros_like(z)
    out[visited] = floor + (1.0 - floor) * (nz / top) ** gamma
    return np.clip(out, 0.0, 1.0)


def build(lams, hs, lam0, exact=False, xlim=(LAM_LO, LAM_HI)):
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.widgets import CheckButtons, Slider

    ink = LinearSegmentedColormap.from_list("ink", ["#ffffff", "#5a6570"])
    fig = plt.figure(figsize=(12.2, 7.4))
    manager = getattr(fig.canvas, "manager", None)
    if manager is not None:
        manager.set_window_title("Core entropy, bifurcation diagram and tent map")

    ax_h = fig.add_axes([0.045, 0.175, 0.415, 0.745])
    ax_b = ax_h.twinx()
    ax_s = fig.add_axes([0.045, 0.085, 0.415, 0.026])
    # The three right-hand panels are placed so that the commuting square can be drawn
    # with straight segments: A directly above g (so a vertical line is the same x), and
    # the tent map to A's right with its argument on the vertical axis (so a horizontal
    # line at height A(x) meets the graph at f_s(A(x)) = A(g(x))).
    ax_g = fig.add_axes([0.525, 0.165, 0.200, 0.270])
    ax_A = fig.add_axes([0.525, 0.520, 0.200, 0.270])
    ax_t = fig.add_axes([0.775, 0.520, 0.200, 0.270])
    ax_cb = fig.add_axes([0.775, 0.395, 0.130, 0.060])
    ax_n = fig.add_axes([0.775, 0.315, 0.130, 0.025])
    ax_w = fig.add_axes([0.775, 0.245, 0.130, 0.025])

    # ---- the bifurcation diagram, behind the entropy curve
    lo, hi = xlim
    # sized to about twice the panel's pixels: more is invisible, and the image is
    # redrawn on every slider event, so it is not free.
    nlam, nrows = 800, 420
    z = bifurcation(lo, hi, nlam=nlam, nrows=nrows)
    im = ax_b.imshow(z, origin="lower", extent=[lo, hi, 0.0, 1.0], aspect="auto",
                     cmap=ink, vmin=0.0, vmax=max(1.0, float(np.percentile(z, 99.5))),
                     interpolation="nearest", zorder=1)
    ax_b.set_ylim(0.0, 1.0)
    ax_b.set_yticks([])
    ax_b.set_xlim(lo, hi)
    ax_b.set_ylabel("attracting set", fontsize=8, color="0.45")
    ax_b.tick_params(labelsize=8)
    ax_h.set_zorder(ax_b.get_zorder() + 1)
    ax_h.patch.set_visible(False)

    pin = {"busy": False}

    def pin_y(_ax=None):
        """Keep the entropy axis where it is: only lambda is zoomed or panned."""
        if pin["busy"]:
            return
        pin["busy"] = True
        try:
            if tuple(ax_h.get_ylim()) != Y_PIN:
                ax_h.set_ylim(*Y_PIN)
        finally:
            pin["busy"] = False

    # ---- the entropy curve, a marker and a vertical line at lambda
    ax_h.plot(lams, hs, color=BLUE, lw=1.3, zorder=5,
              path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])
    for x in (ce.LAM_INF, ce.LAM_P3):
        ax_h.axvline(x, color="0.35", lw=0.8, ls="--", zorder=6)
    ax_h.axhline(LOG2, color="0.45", lw=0.8, ls=":", zorder=6)
    ax_h.set_xlim(lo, hi)
    ax_h.set_ylim(-0.02, 0.70)
    ax_h.set_xlabel(r"$\lambda$", fontsize=10)
    ax_h.set_ylabel(r"$h(\lambda)$", fontsize=10)
    ax_h.set_title(r"core entropy and the bifurcation diagram", fontsize=10)
    ax_h.text(ce.LAM_INF - 0.02, 0.655, r"$\lambda_\infty$", ha="right", va="center",
              fontsize=9, color="0.35", zorder=6, clip_on=True)
    ax_h.text(ce.LAM_P3 + 0.02, 0.655, "period 3", va="center", fontsize=9,
              color="0.35", zorder=6, clip_on=True)
    h0 = ce.h_at(lams, hs, lam0)
    vline = ax_h.axvline(lam0, color=RED, lw=0.9, alpha=0.75, zorder=7)
    mark, = ax_h.plot([lam0], [h0], "o", ms=7, mfc=RED, mec="white", mew=1.0, zorder=8)

    # ---- g_lambda below, and the tent map to the right with its argument on the
    #      vertical axis, so the two panels share a vertical scale with A
    xs = np.linspace(0.0, 1.0, 1201)
    for a in (ax_g, ax_A, ax_t):
        a.set_xlim(0.0, 1.0)
        a.set_ylim(0.0, 1.02)
        a.set_aspect("equal", adjustable="box")
        a.tick_params(labelsize=8)
    ax_g.plot([0.0, 1.0], [0.0, 1.0], color="0.78", lw=0.7, ls=":", zorder=1)
    ax_g.axvline(0.5, color="0.87", lw=0.7, ls="--", zorder=1)
    ax_g.set_xlabel(r"$x$", fontsize=9)
    gline, = ax_g.plot(xs, lam0 * xs * (1.0 - xs), color=BLUE, lw=1.4, zorder=3)
    cw_g, = ax_g.plot([], [], color=ORANGE, lw=0.7, alpha=0.9, zorder=4)
    ax_g.set_ylabel(r"$g_\lambda(x)$", fontsize=9)
    ax_g.set_title(r"$g_\lambda$", fontsize=9, pad=3)

    ax_t.plot([0.0, 1.0], [0.0, 1.0], color="0.78", lw=0.7, ls=":", zorder=1)
    ax_t.axhline(0.5, color="0.87", lw=0.7, ls="--", zorder=1)
    ax_t.set_xlabel(r"$f_s(y)$", fontsize=9)
    ax_t.set_ylabel(r"$y$", fontsize=9)
    ax_t.set_title(r"$f_s$, argument up", fontsize=9, pad=3)
    tline, = ax_t.plot(tent(xs, np.exp(h0)), xs, color=GREEN, lw=1.4, zorder=3)
    cw_t, = ax_t.plot([], [], color=ORANGE, lw=0.7, alpha=0.9, zorder=4)

    # ---- the semi-conjugacy A
    xg = np.linspace(0.0, 1.0, 241)
    aline, = ax_A.plot([], [], color="#7a4b16", lw=1.4, zorder=3)
    adots, = ax_A.plot([], [], "o", ms=3.2, mfc=ORANGE, mec="white", mew=0.4, zorder=4)
    ax_A.plot([0.0, 1.0], [0.0, 1.0], color="0.85", lw=0.7, ls=":", zorder=1)
    ax_A.set_ylabel(r"$A(x)$", fontsize=9)
    ax_A.set_title(r"$A$", fontsize=9, pad=3)
    ax_A.tick_params(labelbottom=False)   # g carries the shared x axis below it

    # ---- the commuting square: up from g, across to the tent map
    # For any x, A(g(x)) = f_s(A(x)).  The vertical segment is a single x read off two
    # graphs; the horizontal one is the value A(x) fed to f_s, and where it meets the
    # tent graph -- drawn with its argument up -- is f_s(A(x)), the same number as A at
    # the image point.  In the cobweb view this is shown at a few points of the orbit,
    # which is what carries the g-orbit onto the f_s-orbit; one per point would be a mess.
    from matplotlib.patches import ConnectionPatch

    X0 = 0.3            # the point used when the plain graphs are on show
    NSQ = 4             # how many orbit points get a connection in the cobweb view
    squares = []
    for k in range(NSQ):
        # a red family, distinct from the orange cobweb it is drawn over
        col = (0.76, 0.05 + 0.10 * k, 0.16, 1.0 - 0.18 * k)
        up = ConnectionPatch((X0, 0.0), (X0, 0.0), coordsA="data", axesA=ax_g,
                             coordsB="data", axesB=ax_A, color=col, lw=0.9, zorder=6)
        across = ConnectionPatch((0.0, 0.0), (0.0, 0.0), coordsA="data", axesA=ax_A,
                                 coordsB="data", axesB=ax_t, color=col, lw=0.9, zorder=6)
        fig.add_artist(up)
        fig.add_artist(across)
        pg, = ax_g.plot([X0], [0.0], "o", ms=3.6, mfc=col, mec="white", mew=0.5, zorder=7)
        pa, = ax_A.plot([X0], [0.0], "o", ms=3.6, mfc=col, mec="white", mew=0.5, zorder=7)
        pt, = ax_t.plot([0.0], [0.0], "o", ms=3.6, mfc=col, mec="white", mew=0.5, zorder=7)
        squares.append({"up": up, "across": across, "arts": (up, across, pg, pa, pt),
                        "pg": pg, "pa": pa, "pt": pt})

    def hide_squares():
        for sq in squares:
            for art in sq["arts"]:
                art.set_visible(False)

    def place_square(k, x0, lam):
        """Put one connection at x0: up from g, across from A to the tent graph."""
        s, sq = cur["s"], squares[k]
        gx = lam * x0 * (1.0 - x0)
        ax_ = float(semiconjugacy([x0], lam, s)[0])
        fx = float(tent(ax_, s))
        sq["pg"].set_data([x0], [gx])
        sq["pa"].set_data([x0], [ax_])
        sq["pt"].set_data([fx], [ax_])
        sq["up"].xy1 = (x0, gx)
        sq["up"].xy2 = (x0, ax_)
        sq["across"].xy1 = (x0, ax_)
        sq["across"].xy2 = (fx, ax_)
        for art in sq["arts"]:
            art.set_visible(True)
        return gx, ax_, fx

    # ---- the widgets
    sl = Slider(ax_s, r"$\lambda$", LAM_LO, LAM_HI, valinit=lam0, valstep=0.0005,
                valfmt="%.4f", color=BLUE)
    sl.label.set_fontsize(10)
    sl.valtext.set_fontsize(9)
    # Slider.set_val draws before it notifies its observers, so leaving its own
    # drawing on would render the figure twice per drag, the first time stale.
    sl.drawon = False

    cb = CheckButtons(ax_cb, ["iterates", "A"], [True, False])
    ax_cb.set_title("view", fontsize=8)
    for lab in cb.labels:
        lab.set_fontsize(8)
    ax_A.set_visible(cb.get_status()[1])   # the checkbox's opening state is the panel's

    sl_w = Slider(ax_w, "cobweb steps", 5, 200, valinit=COBWEB_STEPS, valstep=1,
                  valfmt="%d", color=ORANGE, handle_style={"size": 0.86})
    sl_w.label.set_fontsize(8)
    sl_w.valtext.set_fontsize(8)
    sl_w.drawon = False

    sl_n = Slider(ax_n, r"$n$  (+/-)", 1, 10, valinit=1, valstep=1, valfmt="%d",
                  color=GREEN, handle_style={"size": 0.86})
    sl_n.label.set_fontsize(9)
    sl_n.valtext.set_fontsize(9)
    sl_n.drawon = False

    info = fig.text(0.045, 0.010, "", fontsize=9.0, va="bottom", family="monospace")

    cache = {}
    cur = {"lam": lam0, "h": h0, "s": float(np.exp(h0))}
    _orbit = [0.5] * (NSQ + 1)
    state = {"range": (lo, hi)}
    busy = {"on": False}

    def h_of(lam):
        """h at the slider's lambda, plus how long it took and whether it was cached."""
        if not exact:
            i = int(np.argmin(np.abs(lams - lam)))
            sl.valtext.set_text("%.5f" % float(lams[i]))
            return float(lams[i]), float(hs[i]), None
        key = round(lam, 7)
        if key in cache:
            return lam, cache[key], "cached"
        t0 = time.perf_counter()
        h = float(ce.core_entropy_one(lam))
        cache[key] = h
        return lam, h, "%.0f ms" % (1e3 * (time.perf_counter() - t0))

    def draw_maps():
        """The two map panels: cobweb of the critical orbit, or the n-th iterates."""
        lam, s = cur["lam"], cur["s"]
        n = int(round(sl_n.val))
        if cb.get_status()[0]:
            gl = lam * xs * (1.0 - xs)
            for _ in range(n - 1):
                gl = lam * gl * (1.0 - gl)
            tl = tent(xs, s)
            for _ in range(n - 1):
                tl = tent(tl, s)
            gline.set_ydata(gl)
            tline.set_xdata(tl)          # the tent panel has its argument on y
            cw_g.set_data([], [])
            cw_t.set_data([], [])
            lg = lap_count(lam, s, n, "g")
            lt = lap_count(lam, s, n, "t")
            ax_g.set_title(r"$g_\lambda^{%d}$, $\ell_{%d}=%d$" % (n, n, lg), fontsize=9)
            ax_t.set_title(r"$f_s^{%d}$, $\ell_{%d}=%d$" % (n, n, lt), fontsize=9)
            ax_g.set_ylabel(r"$g_\lambda^{%d}(x)$" % n, fontsize=9)
            ax_t.set_ylabel(r"$f_s^{%d}(y)$" % n, fontsize=9)
        else:
            gline.set_ydata(lam * xs * (1.0 - xs))
            tline.set_xdata(tent(xs, s))     # argument on y, value across
            steps = int(round(sl_w.val))
            cx, cy = cobweb(lambda y: lam * y * (1.0 - y), 0.5, steps)
            cw_g.set_data(cx, cy)
            tx, ty = cobweb(lambda y: tent(y, s), 0.5, steps)
            cw_t.set_data(ty, tx)            # transposed with the panel
            ax_g.set_title(r"$g_\lambda$", fontsize=9)
            ax_t.set_title(r"$f_s$, argument up", fontsize=9)
            ax_g.set_ylabel(r"$g_\lambda(x)$", fontsize=9)
            ax_t.set_ylabel(r"$f_s(y)$", fontsize=9)

    def update(_=None):
        lam, h, timing = h_of(float(sl.val))
        s = float(np.exp(h))
        cur.update(lam=lam, h=h, s=s)
        c = lam / 2.0 - lam ** 2 / 4.0
        _orbit[:] = critical_orbit(lam, NSQ + 1)
        mark.set_data([lam], [h])
        vline.set_xdata([lam, lam])
        draw_maps()
        # The connections are drawn only when the panels really show g and f_s: in the
        # cobweb view at a few points of the orbit, in the iterates view at n = 1 only.
        hide_squares()
        cobweb_view = not cb.get_status()[0]
        plain = cobweb_view or int(round(sl_n.val)) == 1
        shown = []
        if plain and ax_A.get_visible():
            if cobweb_view:
                steps = int(round(sl_w.val))
                for k in range(min(NSQ, steps - 1)):
                    shown.append(place_square(k, 0.5 if k == 0 else _orbit[k], lam))
            else:
                shown.append(place_square(0, X0, lam))
        if ax_A.get_visible():
            o = critical_orbit(lam)
            a = semiconjugacy(xg, lam, s)
            aline.set_data(xg, a)
            adots.set_data(o, semiconjugacy(o, lam, s))
            err = np.abs(semiconjugacy(lam * xg * (1.0 - xg), lam, s) - tent(a, s)).max()
            if shown:
                gx0, ax0, fx0 = shown[0]
                second = (f"A(x) = {ax0:.6f}  f_s(A(x)) = {fx0:.6f}  "
                          f"residual {err:.1e} over {xg.size} points")
            else:
                second = f"residual {err:.1e} over {xg.size} points"
        else:
            second = "A hidden"
        info.set_text(f"lambda = {lam:.5f}    c = {c:.6f}      "
                      f"h = {h:.6f}  h/log2 = {h / LOG2:.4f}  s = {s:.6f}"
                      + (f"    [{timing}]" if exact else "")
                      + "\n" + second)
        fig.canvas.draw_idle()

    def set_slider_range(lo_, hi_):
        """Rescale the slider to the visible lambda window, so its ends are the axis ends."""
        sl.valmin, sl.valmax = lo_, hi_
        sl.valstep = (hi_ - lo_) / 2000.0
        sl.ax.set_xlim(lo_, hi_)
        sl.set_val(min(max(float(sl.val), lo_), hi_))   # clamps and repositions the handle

    def redraw_bif(lo_, hi_, force=False):
        lo_, hi_ = max(LAM_LO, lo_), min(LAM_HI, hi_)
        if hi_ - lo_ < 2e-4:
            return
        plo, phi = state["range"]
        # Skip only when the window still lies inside the computed one and covers at
        # least half of it: a scroll-wheel zoom arrives in many small steps and never
        # fires a button release, so a width-only test would leave a stale, coarse
        # image stretched over the new window.
        if not force and lo_ >= plo and hi_ <= phi and (hi_ - lo_) > 0.55 * (phi - plo):
            return
        z = bifurcation(lo_, hi_, nlam=nlam, nrows=nrows)
        im.set_data(z)
        im.set_extent([lo_, hi_, 0.0, 1.0])
        im.set_clim(0.0, max(1.0, float(np.percentile(z, 99.5))))
        state["range"] = (lo_, hi_)
        set_slider_range(lo_, hi_)
        fig.canvas.draw_idle()

    def on_xlim(_ax):
        if busy["on"]:
            return
        busy["on"] = True
        try:
            a, b = sorted(_ax.get_xlim())
            redraw_bif(a, b)
        finally:
            busy["on"] = False

    ax_h.callbacks.connect("xlim_changed", on_xlim)
    ax_h.callbacks.connect("ylim_changed", pin_y)

    # Zoom and pan act on lambda alone.  The toolbar applies both axes, but it already
    # knows how to apply one: release_zoom passes the held key to _set_view_from_bbox as
    # its mode, and _get_pan_points zeroes the vertical delta for the key "x".  So the
    # two methods the toolbar ends up calling get x-only versions.
    #
    # Both axes need it.  The toolbar collects the axes under the cursor and then adds
    # their twinned siblings unconditionally, and the entropy axes' patch is invisible so
    # that the diagram shows through it -- which means the bifurcation axes is picked up
    # as well, and without this it is the diagram's own y that zooms.
    def x_only(ax):
        bbox, pan = ax._set_view_from_bbox, ax.drag_pan
        ax._set_view_from_bbox = (
            lambda b, direction="in", mode=None, twinx=False, twiny=False:
            bbox(b, direction=direction, mode="x", twinx=twinx, twiny=twiny))
        ax.drag_pan = lambda button, key, x, y: pan(button, "x", x, y)

    for _ax in (ax_h, ax_b):
        x_only(_ax)

    pin_b = {"busy": False}

    def pin_y_b(_ax=None):
        """The attractor axis is the interval [0, 1] and never moves either."""
        if pin_b["busy"]:
            return
        pin_b["busy"] = True
        try:
            if tuple(ax_b.get_ylim()) != (0.0, 1.0):
                ax_b.set_ylim(0.0, 1.0)
        finally:
            pin_b["busy"] = False

    ax_b.callbacks.connect("ylim_changed", pin_y_b)
    fig.canvas.mpl_connect(
        "button_release_event",
        lambda _e: redraw_bif(*sorted(ax_h.get_xlim()), force=True))
    def step_n(delta):
        """Move n by one, clamped to 1..10; the slider follows."""
        sl_n.set_val(max(1, min(10, int(round(sl_n.val)) + delta)))

    def on_key(event):
        if event.key in ("+", "=", "up"):
            step_n(+1)
        elif event.key in ("-", "_", "down"):
            step_n(-1)

    def on_check(_label):
        ax_A.set_visible(cb.get_status()[1])
        update()

    cb.on_clicked(on_check)
    fig.canvas.mpl_connect("key_press_event", on_key)
    # a change of n changes whether the commuting square applies, so update() rather
    # than just draw_maps()
    sl_n.on_changed(lambda _v: update() if cb.get_status()[0] else None)
    sl_w.on_changed(lambda _v: update() if not cb.get_status()[0] else None)

    sl.on_changed(update)
    set_slider_range(lo, hi)
    update()
    return fig, {"slider": sl, "iterates": cb, "n": sl_n, "steps": sl_w, "info": info,
                 "axes": {"h": ax_h, "b": ax_b, "g": ax_g, "t": ax_t, "A": ax_A},
                 "lines": {"g": gline, "t": tline, "mark": mark, "vline": vline,
                           "A": aline, "Adots": adots},
                 "labels": {"g_ylabel": ax_g.yaxis.label, "t_ylabel": ax_t.yaxis.label}}


def start_with_zoom(fig):
    """Select the toolbar's zoom tool, if this backend has a toolbar.

    The toolbar is created with the window, which for some backends happens only at the
    first draw, so this asks once immediately and once more on that draw.
    """
    state = {"done": False}

    def once(_event=None):
        if state["done"]:
            return
        tb = getattr(fig.canvas, "toolbar", None)
        if tb is None or not hasattr(tb, "zoom"):
            return
        if not state.get("band") and hasattr(tb, "draw_rubberband"):
            # the zoom rectangle would otherwise be drawn as a box while only its width
            # is used; stretch it to the height of the axes so the band says what it does
            orig = tb.draw_rubberband
            box = fig.axes[0].bbox

            def band(event, x0, y0, x1, y1, *args, **kwargs):
                return orig(event, x0, box.y0, x1, box.y1, *args, **kwargs)

            tb.draw_rubberband = band
            state["band"] = True
        tb.zoom()
        state["done"] = True

    fig.canvas.mpl_connect("draw_event", once)
    once()


def main():
    ap = argparse.ArgumentParser(
        description="core entropy, bifurcation diagram and the tent map, in a window")
    ap.add_argument("--lam", type=float, default=ce.LAM_P3,
                    help="initial lambda (default: the period-3 centre)")
    ap.add_argument("--exact", action="store_true",
                    help="compute h at the slider's lambda (about 10-36 ms each) "
                         "instead of reading it off the grid")
    ap.add_argument("--xlim", type=float, nargs=2, metavar=("LO", "HI"),
                    help="start with the lambda axis on this range")
    ap.add_argument("--save", metavar="PNG",
                    help="render one frame to this file instead of opening a window")
    args = ap.parse_args()

    import matplotlib
    if args.save:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xlim = tuple(args.xlim) if args.xlim else (LAM_LO, LAM_HI)
    lams, hs = sweep()
    fig, _ = build(lams, hs, args.lam, exact=args.exact, xlim=xlim)

    if args.save:
        fig.savefig(args.save, dpi=130, bbox_inches="tight")
        print(f"wrote {args.save}")
        return
    print(f"h({args.lam:.4f}) = {ce.h_at(lams, hs, args.lam):.6f}; "
          f"drag the slider, zoom the lambda axis")
    start_with_zoom(fig)
    plt.show()


if __name__ == "__main__":
    main()

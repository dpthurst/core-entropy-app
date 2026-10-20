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
points on its graph and the residual max|f_s(A(x)) - A(g(x))| in the readout.  Two
checkboxes switch the iterates view and the A panel; with A off a slider event costs
48 ms rather than 60.  A red vertical line marks the current lambda.

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
BIF_GAMMA, BIF_FLOOR = 0.5, 0.5   # ink = floor + (1-floor)*(count/top)**gamma
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
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.widgets import CheckButtons, Slider

    ink = LinearSegmentedColormap.from_list("ink", ["#ffffff", "#97a1ab"])
    fig = plt.figure(figsize=(12.2, 7.4))

    ax_h = fig.add_axes([0.045, 0.175, 0.415, 0.745])
    ax_b = ax_h.twinx()
    ax_s = fig.add_axes([0.045, 0.085, 0.415, 0.026])
    ax_g = fig.add_axes([0.525, 0.600, 0.200, 0.290])
    ax_t = fig.add_axes([0.775, 0.600, 0.200, 0.290])
    ax_A = fig.add_axes([0.525, 0.180, 0.200, 0.290])
    ax_cb = fig.add_axes([0.775, 0.395, 0.130, 0.055])
    ax_n = fig.add_axes([0.775, 0.300, 0.130, 0.025])

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
    ax_h.plot(lams, hs, color=BLUE, lw=1.3, zorder=5)
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

    # ---- g_lambda and the tent map, with cobwebs or iterates
    xs = np.linspace(0.0, 1.0, 1201)
    for a in (ax_g, ax_t):
        a.plot([0.0, 1.0], [0.0, 1.0], color="0.78", lw=0.7, ls=":", zorder=1)
        a.axvline(0.5, color="0.87", lw=0.7, ls="--", zorder=1)
        a.set_xlim(0.0, 1.0)
        a.set_ylim(0.0, 1.02)
        a.set_aspect("equal", adjustable="box")
        a.set_xlabel(r"$x$", fontsize=9)
        a.tick_params(labelsize=8)
    gline, = ax_g.plot(xs, lam0 * xs * (1.0 - xs), color=BLUE, lw=1.4, zorder=3)
    tline, = ax_t.plot(xs, tent(xs, np.exp(h0)), color=GREEN, lw=1.4, zorder=3)
    cw_g, = ax_g.plot([], [], color=ORANGE, lw=0.7, alpha=0.9, zorder=4)
    cw_t, = ax_t.plot([], [], color=ORANGE, lw=0.7, alpha=0.9, zorder=4)
    ax_g.set_ylabel(r"$g_\lambda(x)$", fontsize=9)
    ax_t.set_ylabel(r"$f_s(y)$", fontsize=9)

    # ---- the semi-conjugacy A
    xg = np.linspace(0.0, 1.0, 241)
    aline, = ax_A.plot([], [], color="#7a4b16", lw=1.4, zorder=3)
    adots, = ax_A.plot([], [], "o", ms=3.2, mfc=ORANGE, mec="white", mew=0.4, zorder=4)
    ax_A.plot([0.0, 1.0], [0.0, 1.0], color="0.85", lw=0.7, ls=":", zorder=1)
    ax_A.set_xlim(0.0, 1.0)
    ax_A.set_ylim(0.0, 1.02)
    ax_A.set_aspect("equal", adjustable="box")
    ax_A.set_xlabel(r"$x$", fontsize=9)
    ax_A.set_ylabel(r"$A(x)$", fontsize=9)
    ax_A.tick_params(labelsize=8)
    ax_A.set_title(r"$A$: $f_s\circ A=A\circ g_\lambda$", fontsize=9)

    # ---- the widgets
    sl = Slider(ax_s, r"$\lambda$", LAM_LO, LAM_HI, valinit=lam0, valstep=0.0005,
                valfmt="%.4f", color=BLUE)
    sl.label.set_fontsize(10)
    sl.valtext.set_fontsize(9)
    # Slider.set_val draws before it notifies its observers, so leaving its own
    # drawing on would render the figure twice per drag, the first time stale.
    sl.drawon = False

    cb = CheckButtons(ax_cb, ["iterates", "A"], [False, True])
    ax_cb.set_title("view", fontsize=8)
    for lab in cb.labels:
        lab.set_fontsize(8)

    sl_n = Slider(ax_n, r"$n$  (+/-)", 1, 10, valinit=4, valstep=1, valfmt="%d",
                  color=GREEN, handle_style={"size": 0.86})
    sl_n.label.set_fontsize(9)
    sl_n.valtext.set_fontsize(9)
    sl_n.drawon = False

    info = fig.text(0.045, 0.010, "", fontsize=9.0, va="bottom", family="monospace")

    cache = {}
    cur = {"lam": lam0, "h": h0, "s": float(np.exp(h0))}
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
            tline.set_ydata(tl)
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
            tline.set_ydata(tent(xs, s))
            cx, cy = cobweb(lambda y: lam * y * (1.0 - y), 0.5, 60)
            cw_g.set_data(cx, cy)
            tx, ty = cobweb(lambda y: tent(y, s), 0.5, 60)
            cw_t.set_data(tx, ty)
            ax_g.set_title(r"$g_\lambda$: critical orbit", fontsize=9)
            ax_t.set_title(r"$f_s$: the same orbit", fontsize=9)
            ax_g.set_ylabel(r"$g_\lambda(x)$", fontsize=9)
            ax_t.set_ylabel(r"$f_s(y)$", fontsize=9)

    def update(_=None):
        lam, h, timing = h_of(float(sl.val))
        s = float(np.exp(h))
        cur.update(lam=lam, h=h, s=s)
        c = lam / 2.0 - lam ** 2 / 4.0
        mark.set_data([lam], [h])
        vline.set_xdata([lam, lam])
        draw_maps()
        if ax_A.get_visible():
            o = critical_orbit(lam)
            a = semiconjugacy(xg, lam, s)
            aline.set_data(xg, a)
            adots.set_data(o, semiconjugacy(o, lam, s))
            err = np.abs(semiconjugacy(lam * xg * (1.0 - xg), lam, s) - tent(a, s)).max()
            second = f"max |f_s(A(x)) - A(g_lambda(x))| = {err:.1e}  over {xg.size} points"
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
    _bbox, _pan = ax_h._set_view_from_bbox, ax_h.drag_pan
    ax_h._set_view_from_bbox = (
        lambda bbox, direction="in", mode=None, twinx=False, twiny=False:
        _bbox(bbox, direction=direction, mode="x", twinx=twinx, twiny=twiny))
    ax_h.drag_pan = lambda button, key, x, y: _pan(button, "x", x, y)
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
    sl_n.on_changed(lambda _v: (draw_maps(), fig.canvas.draw_idle())
                    if cb.get_status()[0] else None)

    sl.on_changed(update)
    set_slider_range(lo, hi)
    update()
    return fig, {"slider": sl, "iterates": cb, "n": sl_n, "info": info,
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

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
scales readable.  Right: g_lambda and the tent map f_s of slope s = exp(h(lambda)),
each with the cobweb of the critical orbit.

By default h is read off the 28 696-parameter grid, which is instant; with --exact
it is computed at the slider's own lambda by core_entropy.py, at about 10 ms above
lambda = 3.6 and 36 ms below, where the solver must quadruple N before it can
report zero.  A full redraw is a further 25 ms.  Results are cached, so revisiting
a value is free.  A recomputed bifurcation diagram costs about 10-15 ms, and is
built once per zoom rather than per mouse movement.

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


def bifurcation(lo, hi, nlam=640, nrows=400, transient=1200, keep=700):
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


def build(lams, hs, lam0, exact=False, xlim=(LAM_LO, LAM_HI)):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.widgets import Slider

    ink = LinearSegmentedColormap.from_list("ink", ["#ffffff", "#97a1ab"])
    fig = plt.figure(figsize=(11.5, 6.2))

    ax_h = fig.add_axes([0.055, 0.13, 0.49, 0.78])
    ax_b = ax_h.twinx()
    ax_g = fig.add_axes([0.60, 0.59, 0.35, 0.34])
    ax_t = fig.add_axes([0.60, 0.21, 0.35, 0.34])
    ax_s = fig.add_axes([0.60, 0.12, 0.35, 0.03])

    # ---- the bifurcation diagram, behind the entropy curve
    lo, hi = xlim
    z = bifurcation(lo, hi)
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

    # ---- the entropy curve
    ax_h.plot(lams, hs, color=BLUE, lw=1.3, zorder=5)
    for x in (ce.LAM_INF, ce.LAM_P3):
        ax_h.axvline(x, color="0.35", lw=0.8, ls="--", zorder=6)
    ax_h.axhline(LOG2, color="0.45", lw=0.8, ls=":", zorder=6)
    ax_h.set_xlim(lo, hi)
    ax_h.set_ylim(-0.02, 0.70)
    ax_h.set_xlabel(r"$\lambda$", fontsize=10)
    ax_h.set_ylabel(r"$h(\lambda)$", fontsize=10)
    ax_h.set_title(r"core entropy and the bifurcation diagram,   "
                   r"$g_\lambda(x)=\lambda x(1-x)$", fontsize=10)
    ax_h.text(ce.LAM_INF - 0.02, 0.655, r"$\lambda_\infty$", ha="right", va="center",
              fontsize=9, color="0.35", zorder=6, clip_on=True)
    ax_h.text(ce.LAM_P3 + 0.02, 0.655, "period 3", va="center", fontsize=9,
              color="0.35", zorder=6, clip_on=True)
    h0 = ce.h_at(lams, hs, lam0)
    mark, = ax_h.plot([lam0], [h0], "o", ms=7, mfc=RED, mec="white", mew=1.0, zorder=8)

    # ---- the two maps, each with the cobweb of the critical orbit
    xs = np.linspace(0.0, 1.0, 601)
    for a in (ax_g, ax_t):
        a.plot([0.0, 1.0], [0.0, 1.0], color="0.78", lw=0.7, ls=":", zorder=1)
        a.axvline(0.5, color="0.87", lw=0.7, ls="--", zorder=1)
        a.set_xlim(0.0, 1.0)
        a.set_ylim(0.0, 1.02)
        a.set_aspect("equal", adjustable="box")
        a.set_xlabel(r"$x$", fontsize=9)
        a.tick_params(labelsize=8)
    gline, = ax_g.plot(xs, lam0 * xs * (1.0 - xs), color=BLUE, lw=1.5, zorder=3,
                       label=r"$g_\lambda$")
    tline, = ax_t.plot(xs, tent(xs, np.exp(h0)), color=GREEN, lw=1.5, zorder=3,
                       label=r"$f_s,\ s=e^{h(\lambda)}$")
    cw_g, = ax_g.plot([], [], color=ORANGE, lw=0.7, alpha=0.9, zorder=4,
                      label="critical orbit")
    cw_t, = ax_t.plot([], [], color=ORANGE, lw=0.7, alpha=0.9, zorder=4,
                      label="critical orbit")
    for a in (ax_g, ax_t):
        a.legend(loc="upper left", fontsize=7.5, frameon=False)
    ax_g.set_ylabel(r"$g_\lambda(x)$", fontsize=9)
    ax_t.set_ylabel(r"$f_s(y)$", fontsize=9)
    ax_g.set_title(r"$g_\lambda$ and the cobweb of the critical orbit", fontsize=9)
    ax_t.set_title(r"the tent map $f_s$, with its own cobweb", fontsize=9)

    info = fig.text(0.60, 0.035, "", fontsize=9.5, va="bottom", family="monospace")

    sl = Slider(ax_s, r"$\lambda$", LAM_LO, LAM_HI, valinit=lam0, valstep=0.0005,
                valfmt="%.4f", color=BLUE)
    sl.label.set_fontsize(10)
    sl.valtext.set_fontsize(9)

    cache = {}
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

    def update(_=None):
        lam, h, timing = h_of(float(sl.val))
        s = float(np.exp(h))
        c = lam / 2.0 - lam ** 2 / 4.0
        mark.set_data([lam], [h])
        gline.set_ydata(lam * xs * (1.0 - xs))
        tline.set_ydata(tent(xs, s))
        cx, cy = cobweb(lambda y: lam * y * (1.0 - y), 0.5, 60)
        cw_g.set_data(cx, cy)
        tx, ty = cobweb(lambda y: tent(y, s), 0.5, 60)
        cw_t.set_data(tx, ty)
        info.set_text(f"lambda = {lam:.5f}    c = {c:.6f}\n"
                      f"h = {h:.6f}    h/log2 = {h / LOG2:.4f}    s = {s:.6f}"
                      + (f"    [{timing}]" if exact else ""))
        fig.canvas.draw_idle()

    def redraw_bif(lo_, hi_, force=False):
        lo_, hi_ = max(LAM_LO, lo_), min(LAM_HI, hi_)
        if hi_ - lo_ < 2e-4:
            return
        plo, phi = state["range"]
        if not force and abs(np.log((hi_ - lo_) / (phi - plo))) < 0.06:
            return
        z = bifurcation(lo_, hi_)
        im.set_data(z)
        im.set_extent([lo_, hi_, 0.0, 1.0])
        im.set_clim(0.0, max(1.0, float(np.percentile(z, 99.5))))
        state["range"] = (lo_, hi_)
        m = (lams >= lo_) & (lams <= hi_)
        if m.any():
            top = float(hs[m].max())
            ax_h.set_ylim(-0.02, max(top * 1.15, 0.02) if hi_ - lo_ < 1.9 else 0.70)

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
    fig.canvas.mpl_connect(
        "button_release_event",
        lambda _e: redraw_bif(*sorted(ax_h.get_xlim()), force=True))

    sl.on_changed(update)
    update()
    return fig, sl


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
    plt.show()


if __name__ == "__main__":
    main()

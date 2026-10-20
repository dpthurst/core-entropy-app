#!/usr/bin/env python3
"""
tent_map_app.py -- the core entropy of the logistic family in a native window.

The matplotlib counterpart of tent_map.html (and of core_entropy_makie.jl on
the Julia side): drag the slider to move lambda and both panels follow.  No
browser, no plotly, no Julia -- the only requirement beyond numpy is
matplotlib.

    python3 tent_map_app.py                            # open the window
    python3 tent_map_app.py --lam 3.7                  # start elsewhere
    python3 tent_map_app.py --save frame.png --lam 3.7 # render one frame, headless

Left: h(lambda) from core_entropy.py, with a marker at the current lambda, the
dashed lines at the Feigenbaum accumulation and the period-3 centre, and the
dotted line at log 2.  Right: g_lambda on [0,1] with its critical orbit, and
the tent map f_s of slope s = exp(h(lambda)), which carries the same entropy.

NOT TESTED INTERACTIVELY: written in an environment with no display, so the
window and the slider have never been driven by a real GUI backend.  The
drawing code is exercised by --save, and the entropy it plots is tested by
`make check`.
"""
import argparse
import os

import numpy as np

import core_entropy as ce

HERE = os.path.dirname(os.path.abspath(__file__))
LOG2 = np.log(2.0)
BLUE, GREEN, RED, ORANGE = "#1c5f9e", "#2f6b4f", "#c1272d", "#e08a1e"


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


def critical_orbit(lam, n=400):
    """The first n points of the orbit of 1/2 under g_lambda."""
    x = np.empty(n)
    x[0] = 0.5
    for j in range(1, n):
        x[j] = lam * x[j - 1] * (1.0 - x[j - 1])
    return x


def build(lams, hs, lam0):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider

    fig = plt.figure(figsize=(11.0, 5.4))
    fig.suptitle(r"core entropy of $g_\lambda(x)=\lambda x(1-x)$", fontsize=12)

    ax_h = fig.add_axes([0.065, 0.24, 0.47, 0.64])
    ax_m = fig.add_axes([0.625, 0.34, 0.31, 0.58])
    ax_s = fig.add_axes([0.625, 0.20, 0.31, 0.035])

    # ---- the entropy curve
    ax_h.plot(lams, hs, color=BLUE, lw=1.4)
    for x in (ce.LAM_INF, ce.LAM_P3):
        ax_h.axvline(x, color="0.6", lw=0.8, ls="--")
    ax_h.axhline(LOG2, color="0.7", lw=0.8, ls=":")
    ax_h.set_xlim(2.0, 4.0)
    ax_h.set_ylim(-0.02, 0.70)
    ax_h.set_xlabel(r"$\lambda$", fontsize=10)
    ax_h.set_ylabel(r"$h(\lambda)$", fontsize=10)
    ax_h.tick_params(labelsize=9)
    ax_h.text(ce.LAM_INF - 0.02, 0.655, r"$\lambda_\infty$", ha="right", va="center",
              fontsize=9, color="0.45")
    ax_h.text(ce.LAM_P3 + 0.02, 0.655, "period 3", va="center", fontsize=9, color="0.45")
    h0 = ce.h_at(lams, hs, lam0)
    mark, = ax_h.plot([lam0], [h0], "o", ms=7, mfc=RED, mec="white", mew=1.0, zorder=5)

    # ---- the map, the tent of the same entropy, and the critical orbit
    xs = np.linspace(0.0, 1.0, 601)
    gline, = ax_m.plot(xs, lam0 * xs * (1.0 - xs), color=BLUE, lw=1.8,
                       label=r"$g_\lambda$")
    tline, = ax_m.plot(xs, tent(xs, np.exp(h0)), color=GREEN, lw=1.4,
                       label=r"$f_s,\ s=e^{h(\lambda)}$")
    ax_m.plot([0.0, 1.0], [0.0, 1.0], color="0.7", lw=0.8, ls=":")
    ax_m.axvline(0.5, color="0.8", lw=0.8, ls="--")
    o = critical_orbit(lam0)
    orb, = ax_m.plot(o[:-1], o[1:], "o", ms=3.0, color=ORANGE, alpha=0.7,
                     label="critical orbit")
    ax_m.set_xlim(0.0, 1.0)
    ax_m.set_ylim(0.0, 1.02)
    ax_m.set_aspect("equal", adjustable="box")
    ax_m.set_xlabel(r"$x$", fontsize=10)
    ax_m.set_ylabel(r"$g_\lambda(x)$", fontsize=10)
    ax_m.tick_params(labelsize=9)
    ax_m.legend(loc="upper left", fontsize=8, frameon=False)

    info = fig.text(0.625, 0.06, "", fontsize=9.5, va="bottom", family="monospace")

    sl = Slider(ax_s, r"$\lambda$", 2.0, 4.0, valinit=lam0, valstep=0.0005,
                valfmt="%.4f", color=BLUE)
    sl.label.set_fontsize(10)
    sl.valtext.set_fontsize(9)

    def update(_=None):
        # h is known only at the grid parameters, so snap the slider to one and say so:
        # near the merging parameter the grid spacing alone moves h by 3e-4.
        i = int(np.argmin(np.abs(lams - float(sl.val))))
        lam, h = float(lams[i]), float(hs[i])
        sl.valtext.set_text("%.5f" % lam)
        s = float(np.exp(h))
        c = lam / 2.0 - lam ** 2 / 4.0
        mark.set_data([lam], [h])
        gline.set_ydata(lam * xs * (1.0 - xs))
        tline.set_ydata(tent(xs, s))
        o = critical_orbit(lam)
        orb.set_data(o[:-1], o[1:])
        info.set_text(f"lambda = {lam:.5f}    c = {c:.6f}\n"
                      f"h = {h:.6f}    h/log2 = {h / LOG2:.4f}    s = {s:.6f}")
    # the slider sits on a grid parameter, so what is shown is exactly what was computed
        fig.canvas.draw_idle()

    sl.on_changed(update)
    update()
    return fig, sl


def main():
    ap = argparse.ArgumentParser(
        description="core entropy, critical orbit and the tent map, in a window")
    ap.add_argument("--lam", type=float, default=ce.LAM_P3,
                    help="initial lambda (default: the period-3 centre)")
    ap.add_argument("--save", metavar="PNG",
                    help="render one frame to this file instead of opening a window")
    args = ap.parse_args()

    import matplotlib
    if args.save:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lams, hs = sweep()
    fig, _ = build(lams, hs, args.lam)

    if args.save:
        fig.savefig(args.save, dpi=130, bbox_inches="tight")
        print(f"wrote {args.save}")
        return
    print(f"h({args.lam:.4f}) = {ce.h_at(lams, hs, args.lam):.6f}; drag the slider to move")
    plt.show()


if __name__ == "__main__":
    main()

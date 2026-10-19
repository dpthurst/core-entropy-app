#!/usr/bin/env python3
"""
figs/make_previews.py -- the static previews used in the README.

    python3 figs/make_previews.py

writes preview.png, preview_bifurcation.png, preview_tent.png and
preview_kneading.png in the repository root.  Inputs: core_entropy_data.csv and
figs/fig_kneading.pdf.
"""
import os
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)

plt.rcParams.update({"font.family": "serif",
                     "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
                     "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
                     "font.size": 9})
BLUE, GREEN, RED, GREY = "#1c5f9e", "#2f6b4f", "#c0392b", "#b9c0c7"

d = np.loadtxt("core_entropy_data.csv", delimiter=",", skiprows=1)
LAM, H = d[:, 0], d[:, 2]
L2, LTAU = np.log(2.0), np.log((1 + 5 ** 0.5) / 2)
INF, P3 = 3.5699456718, 3.8318740553
LAM0 = 3.7
S0 = float(np.exp(np.interp(LAM0, LAM, H)))


def g(x, lam=LAM0):
    return lam * x * (1 - x)


def tent(y, s=S0):
    return s * y if y < 0.5 else s * (1 - y)


def cobweb(f, n=34, x0=0.5):
    xs, ys = [], []
    x = x0
    for _ in range(n):
        y = f(x)
        xs += [x, x, y]
        ys += [x, y, y]
        x = y
    return xs, ys


def A_of(x, lam=LAM0, s=S0, n=64):
    """Milnor-Thurston's semiconjugacy: pull a seed back through the branches."""
    y = x
    w = []
    for _ in range(n):
        w.append(1 if y < 0.5 else -1)
        y = g(y, lam)
    u = 0.5
    for c in reversed(w):
        u = u / s if c > 0 else 1 - u / s
    return u


def dedupe(v, tol=1e-9):
    v = sorted(v)
    out = []
    for z in v:
        if not out or z - out[-1] > tol:
            out.append(z)
    return out


def lap_count(pre, n):
    S = [0.5]
    allp = [0.5]
    for _ in range(n - 1):
        nx = []
        for y in S:
            nx += pre(y)
        S = dedupe([p for p in nx if 0 < p < 1])
        allp = dedupe(allp + S)
    return len([p for p in allp if 0 < p < 1]) + 1


def pre_g(y):
    disc = 1 - 4 * y / LAM0
    if disc <= 0:
        return []
    r = np.sqrt(disc)
    return [0.5 * (1 + r), 0.5 * (1 - r)]


def pre_tent(y):
    return [y / S0, 1 - y / S0] if y <= S0 / 2 else []


def density(x0, x1, nc=1000, nr=400, transient=1500, keep=1200):
    lams = np.linspace(x0, x1, nc)
    x = np.full(nc, 0.5)
    for _ in range(transient):
        x = lams * x * (1 - x)
    Z = np.zeros((nr, nc))
    idx = np.arange(nc)
    for _ in range(keep):
        x = lams * x * (1 - x)
        r = np.clip((x * nr).astype(int), 0, nr - 1)
        np.add.at(Z, (r, idx), 1)
    return np.log1p(Z) / np.log1p(Z.max())


INK = mcolors.LinearSegmentedColormap.from_list(
    "ink", [(1, 1, 1, 0), (0.24, 0.24, 0.24, 0.30), (0.10, 0.10, 0.10, 0.82), (0, 0, 0, 0.95)])


def style(a):
    a.tick_params(labelsize=8)
    a.grid(alpha=.2, lw=.5)
    a.spines[["top", "right"]].set_visible(False)


def preview_entropy():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.0, 4.3), dpi=100, layout="constrained",
                                 gridspec_kw={"width_ratios": [1, 1.45], "wspace": 0.16})
    for a in (a1, a2):
        a.plot(LAM, H, lw=1.1, color=BLUE)
        a.set_xlabel(r"$\lambda$", fontsize=11)
        a.set_ylabel(r"core entropy  $h(\lambda)$", fontsize=11)
        a.set_ylim(-0.02, 0.70)
        style(a)
    a1.set_xlim(2, 4)
    for x, lab in ((INF, r"$\lambda_\infty$"), (P3, "period 3")):
        a1.axvline(x, color="0.72", lw=.9, ls=":")
        a1.text(x, 0.67, lab, rotation=90, fontsize=8, color="0.45", ha="right", va="top")
    a1.axhline(L2, color="0.85", lw=.8)
    a1.text(2.03, L2 - 0.045, r"$\log 2$", fontsize=8, color="0.45")
    a1.set_title("full range", fontsize=10.5, color="0.25")
    a2.set_xlim(3.55, 4)
    a2.axhline(L2, color="0.85", lw=.8)
    a2.axhline(LTAU, color="0.85", lw=.8)
    a2.text(3.555, L2 + 0.008, r"$\log 2$", fontsize=8, color="0.45")
    a2.text(3.555, LTAU + 0.008, r"$\log\tau$", fontsize=8, color="0.45")
    a2.axvline(INF, color="0.72", lw=.9, ls=":")
    a2.annotate("period-3 plateau", xy=(3.845, LTAU), xytext=(3.79, 0.60), fontsize=8,
                color="0.35", ha="center",
                arrowprops=dict(arrowstyle="-", color="0.75", lw=.7))
    a2.set_title("the staircase", fontsize=10.5, color="0.25")
    plt.savefig("preview.png", dpi=100, facecolor="white")
    plt.close(fig)


def preview_bifurcation():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.0, 4.4), dpi=100, layout="constrained",
                                 gridspec_kw={"width_ratios": [1, 1.45], "wspace": 0.12})
    for ax, (lo, hi), title in ((a1, (2, 4), "full range"), (a2, (3.72, 4), "the staircase")):
        Z = density(lo, hi)
        t = ax.twinx()
        t.imshow(Z, origin="lower", extent=[lo, hi, 0, 1], aspect="auto", cmap=INK,
                 vmin=0, vmax=1, interpolation="bilinear", zorder=1)
        t.set_ylim(0, 1)
        t.set_yticks([0, 0.5, 1])
        t.tick_params(labelsize=8, colors="0.35")
        t.set_ylabel("x  (limit set)", fontsize=9, color="0.35")
        for sp in t.spines.values():
            sp.set_visible(False)
        ax.plot(LAM, H, lw=1.5, color=BLUE, zorder=5)
        ax.set_zorder(5)
        ax.patch.set_visible(False)
        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.02, 0.70)
        ax.set_xlabel(r"$\lambda$", fontsize=10)
        ax.set_title(title, fontsize=10.5, color="0.25")
        style(ax)
        for xx in (INF, P3):
            if lo < xx < hi:
                ax.axvline(xx, color="0.45", lw=.8, ls=":")
    a1.set_ylabel(r"core entropy  $h(\lambda)$", fontsize=10)
    a2.axhline(L2, color="0.8", lw=.8)
    a2.axhline(LTAU, color="0.8", lw=.8)
    a2.text(3.725, L2 - 0.045, r"$\log 2$", fontsize=8, color="0.4")
    a2.text(3.725, LTAU + 0.008, r"$\log\tau$", fontsize=8, color="0.4")
    plt.savefig("preview_bifurcation.png", dpi=100, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def preview_tent(N=6):
    fig = plt.figure(figsize=(12.6, 6.4), dpi=100)
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.45, wspace=0.28,
                          left=0.06, right=0.975, top=0.94, bottom=0.09)
    ax = fig.add_subplot(gs[0, :])
    Z = density(2, 4)
    t = ax.twinx()
    t.imshow(Z, origin="lower", extent=[2, 4, 0, 1], aspect="auto", cmap=INK, vmin=0, vmax=1,
             interpolation="bilinear", zorder=1)
    t.set_ylim(0, 1)
    t.set_yticks([0, 0.5, 1])
    t.tick_params(labelsize=8, colors="0.35")
    t.set_ylabel("x  (limit set)", fontsize=9, color="0.35")
    for sp in t.spines.values():
        sp.set_visible(False)
    ax.plot(LAM, H, color=BLUE, lw=1.3, zorder=5)
    ax.plot([LAM0], [np.interp(LAM0, LAM, H)], "o", ms=8, color="#c1272d", zorder=6,
            mec="white", mew=1.2)
    ax.set_zorder(5)
    ax.patch.set_visible(False)
    ax.set_xlim(2, 4)
    ax.set_ylim(-0.02, 0.70)
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"core entropy  $h(\lambda)$")
    style(ax)

    x = np.linspace(0, 1, 3000)
    ax = fig.add_subplot(gs[1, 0])
    y = x.copy()
    for _ in range(N):
        y = g(y)
    ax.plot(x, y, color=BLUE, lw=0.7)
    ax.set_title(r"$g_\lambda^{%d}$   with $\ell_{%d}=%d$"
                 % (N, N, lap_count(pre_g, N)), fontsize=9.5)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$g_\lambda^{%d}(x)$" % N)

    ax = fig.add_subplot(gs[1, 1])
    y = x.copy()
    for _ in range(N):
        y = np.where(y < 0.5, S0 * y, S0 * (1 - y))
    ax.plot(x, y, color=GREEN, lw=0.7)
    ax.set_title(r"$f_s^{%d}$   with $\ell_{%d}=%d$"
                 % (N, N, lap_count(pre_tent, N)), fontsize=9.5)
    ax.set_xlabel(r"$y$")
    ax.set_ylabel(r"$f_s^{%d}(y)$" % N)

    ax = fig.add_subplot(gs[1, 2])
    gx = np.linspace(0, 1, 400)
    ax.plot(gx, [A_of(v) for v in gx], color="#333", lw=1.3)
    oL = [0.5]
    for _ in range(6):
        oL.append(g(oL[-1]))
    ax.plot(oL, [A_of(v) for v in oL], "o", ms=5, color=RED, zorder=5)
    ax.plot([0, 1], [0, 1], color="0.75", lw=0.8, ls=":")
    ax.set_title(r"$A$:  $f_s\circ A = A\circ g_\lambda$", fontsize=9.5)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$A(x)$")

    for a in [fig.axes[i] for i in (1, 2, 3)]:
        a.set_xlim(0, 1)
        a.set_ylim(-0.02, 1.04)
        style(a)
    plt.savefig("preview_tent.png", dpi=100, facecolor="white")
    plt.close(fig)


def preview_kneading():
    subprocess.run(["pdftoppm", "-r", "150", "-png", "-singlefile",
                    os.path.join(HERE, "fig_kneading.pdf"), "preview_kneading"],
                   check=True)


if __name__ == "__main__":
    preview_entropy()
    preview_bifurcation()
    preview_tent()
    preview_kneading()
    for f in ("preview.png", "preview_bifurcation.png", "preview_tent.png",
              "preview_kneading.png"):
        print(f"wrote {f} ({os.path.getsize(f)/1024:.0f} KiB)")

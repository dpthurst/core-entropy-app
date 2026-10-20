#!/usr/bin/env python3
"""
Figures for core_entropy_notes.tex.

    python3 figs/make_figures.py        # writes figs/fig*.pdf

The entropy curve is read from core_entropy_data.csv (written by
core_entropy.jl); everything else is recomputed here from the map itself, so
the figures can be rebuilt without Julia.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6, "lines.linewidth": 0.9,
    "savefig.dpi": 600,
})
BLUE, RED, GREY = "#1c5f9e", "#c0392b", "#7f8c8d"
L2, LTAU = np.log(2.0), np.log((1 + 5**0.5) / 2)
LAM_INF, LAM_P3 = 3.5699456718, 3.8318740553

d = np.loadtxt("core_entropy_data.csv", delimiter=",", skiprows=1)
lam, h = d[:, 0], d[:, 2]


def dedupe(v, tol=1e-9):
    """sorted copy with points closer than tol merged"""
    v = sorted(v)
    out = []
    for z in v:
        if not out or z - out[-1] > tol:
            out.append(z)
    return out


def g(x, l):
    return l * x * (1 - x)


def orbit(l, n, x0=0.5):
    xs = [x0]
    for _ in range(n):
        xs.append(g(xs[-1], l))
    return np.array(xs)


def h_of(l):
    return np.interp(l, lam, h)


# ---------------------------------------------------------------- figure 1
def fig_output():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.35),
                                 gridspec_kw={"width_ratios": [1, 1.4], "wspace": 0.22})
    cmap = ink_cmap()
    ylo, yhi = -0.02, 0.74
    # put x = 1 of the density at the same height as h = log 2
    top = 1.0 / ((L2 - ylo) / (yhi - ylo))
    density_axes(a1, 2.0, 4.0, cmap, label=False, top=top)
    density_axes(a2, 3.55, 4.0, cmap, label=True, top=top)
    for a in (a1, a2):
        a.plot(lam, h, color=BLUE, lw=1.0, zorder=5)
        a.set_xlabel(r"$\lambda$")
        a.set_ylim(ylo, yhi)
        a.grid(alpha=0.18, lw=0.5)
        a.spines[["top", "right"]].set_visible(False)
    a1.set_xlim(2, 4)
    a1.set_ylabel(r"core entropy  $h(\lambda)$")
    box = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85)
    a1.axhline(L2, color=GREY, lw=0.6, ls="--")
    a1.text(2.04, L2 + 0.02, r"$\log 2$", color="0.3", fontsize=7, bbox=box)
    for x, t in ((LAM_INF, r"$\lambda_\infty$"), (LAM_P3, r"$3.8319$")):
        a1.axvline(x, color=GREY, lw=0.6, ls=":")
        a1.text(x + 0.015, 0.705, t, ha="left", va="bottom", fontsize=7,
                color="0.25", bbox=box)
    a1.set_title("(a)  the whole range", fontsize=8)
    a2.set_xlim(3.55, 4)
    a2.axhline(L2, color=GREY, lw=0.6, ls="--")
    a2.axhline(LTAU, color=GREY, lw=0.6, ls="--")
    a2.text(3.56, L2 + 0.02, r"$\log 2$", color="0.3", fontsize=7, bbox=box)
    a2.text(3.56, LTAU + 0.02, r"$\log\tau$", color="0.3", fontsize=7, bbox=box)
    for x in (LAM_INF, LAM_P3):
        a2.axvline(x, color=GREY, lw=0.6, ls=":")
    a2.annotate(r"$\lambda_\infty$", xy=(LAM_INF, 0.01), xytext=(3.60, 0.05), fontsize=7,
                color="0.25", bbox=box,
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.5))
    a2.annotate("period-3 plateau", xy=(3.845, LTAU), xytext=(3.80, 0.705), fontsize=7,
                color="0.25", bbox=box,
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.5))
    a2.set_title(r"(b)  the staircase, $\lambda\in[3.55,4]$", fontsize=8)
    fig.savefig("figs/fig_output.pdf", bbox_inches="tight")


# ---------------------------------------------------------------- figure 2
def sign_data(l, n, tol=1e-9):
    """itinerary of the critical value: 0 = left of 1/2, 1 = right, * (a return to
    1/2, within tol) read as 1; and the sign sequence s_k = sign g'(x_k)"""
    xs = orbit(l, n)
    left = (xs[1:] < 0.5) & (np.abs(xs[1:] - 0.5) > tol)
    nu = (~left).astype(int)                 # 1 = right, and a return is read as 1
    s = np.where(left, 1.0, -1.0)
    eps = np.cumprod(np.concatenate([[1.0], s]))
    return xs, nu, eps


def strip(ax, l, n=36):
    xs, nu, eps = sign_data(l, n)
    for k, v in enumerate(nu):
        ax.add_patch(Rectangle((k, 1.15), 0.92, 0.7,
                               color=("#9ecae1" if v == 0 else "#f4a582"), lw=0))
    for k, v in enumerate(eps[1:n + 1]):
        ax.add_patch(Rectangle((k, 0.15), 0.92, 0.7,
                               color=("0.85" if v > 0 else "0.35"), lw=0))
    ax.set_xlim(-0.2, n + 0.2); ax.set_ylim(0, 2.3)
    ax.set_yticks([0.5, 1.5]); ax.set_yticklabels([r"$\varepsilon_k$", r"$\nu_k$"])
    ax.set_xticks(range(0, n + 1, 6))
    ax.tick_params(length=2)
    for sp in ax.spines.values(): sp.set_visible(False)


def fig_kneading():
    lam0 = 3.7
    fig = plt.figure(figsize=(6.6, 4.1))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.45, wspace=0.22)

    # (a) cobweb
    a = fig.add_subplot(gs[0, 0])
    xs = np.linspace(0, 1, 600)
    a.plot(xs, g(xs, lam0), color=BLUE, lw=1.0)
    a.plot([0, 1], [0, 1], color=GREY, lw=0.7, ls=":")
    a.axvline(0.5, color=GREY, lw=0.6, ls=":")
    orb = orbit(lam0, 5)
    x = orb[0]
    for k in range(5):
        y = g(x, lam0)
        a.plot([x, x], [x, y], color=RED, lw=0.7)
        a.plot([x, y], [y, y], color=RED, lw=0.7)
        x = y
    for k, xk in enumerate(orb[:6]):
        a.plot([xk], [g(xk, lam0)], "o", ms=2.6, color=RED)
        a.annotate(str(k), (xk, g(xk, lam0)), textcoords="offset points", xytext=(3, 3),
                   fontsize=6, color=RED)
    a.set_xlim(0, 1); a.set_ylim(0, 1.02)
    a.set_xlabel(r"$x$"); a.set_ylabel(r"$g_\lambda(x)$")
    a.set_title(r"(a)  critical orbit, $\lambda=3.7$", fontsize=8)
    a.text(0.03, 0.93, r"$x_0=\frac{1}{2}$", fontsize=7, color=GREY)

    # (b) the orbit values with their symbols
    a = fig.add_subplot(gs[0, 1])
    m = 24
    xs2, nu2, _ = sign_data(lam0, m)
    k = np.arange(m + 1)
    a.plot(k, xs2, color=RED, lw=0.6, marker="o", ms=2.2, mfc="white", mew=0.7)
    a.axhline(0.5, color=GREY, lw=0.7, ls=":")
    for j in range(1, m + 1):
        a.annotate("L" if xs2[j] < 0.5 else "R", (k[j], xs2[j]), textcoords="offset points",
                   xytext=(0, 5), ha="center", fontsize=5.5, color=BLUE)
    a.set_xlim(-0.5, m + 0.5); a.set_ylim(-0.05, 1.12)
    a.set_xlabel(r"$k$"); a.set_ylabel(r"$x_k$")
    a.set_title(r"(b)  its symbols: the itinerary $\nu$", fontsize=8)
    a.text(0.5, 0.52, r"$x=\frac{1}{2}$", fontsize=6.5, color=GREY, va="bottom")

    # (c),(d) the sequence strips
    a = fig.add_subplot(gs[1, 0]); strip(a, lam0)
    a.set_title(r"(c)  signs for $\lambda=3.7$ (not periodic)", fontsize=8)
    a = fig.add_subplot(gs[1, 1]); strip(a, LAM_P3)
    a.set_title(r"(d)  signs at the period-3 centre", fontsize=8)

    handles = [Rectangle((0, 0), 1, 1, color="#9ecae1"), Rectangle((0, 0), 1, 1, color="#f4a582"),
               Rectangle((0, 0), 1, 1, color="0.85"), Rectangle((0, 0), 1, 1, color="0.35")]
    fig.legend(handles, [r"$\nu_k=0$ (left)", r"$\nu_k=1$ (right)",
                         r"$\varepsilon_k=+1$", r"$\varepsilon_k=-1$"],
               loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.savefig("figs/fig_kneading.pdf", bbox_inches="tight")


# ---------------------------------------------------------------- figure 3
def preimages_log(lam0):
    def f(y):
        disc = 1 - 4*y/lam0
        if disc <= 0:
            return []
        r = np.sqrt(disc)
        return [0.5*(1 + r), 0.5*(1 - r)]
    return f


def preimages_tent(s):
    def f(y):
        if y > s/2:
            return []
        return [p for p in (y/s, 1 - y/s) if 0 < p < 1]
    return f


def lap_count(pre, n, tol=1e-9):
    S = [0.5]; allp = [0.5]
    for _ in range(n - 1):
        nx = []
        for y in S:
            nx += pre(y)
        S = dedupe([p for p in nx if 0 < p < 1], tol)
        allp = dedupe(allp + S, tol)
    return len([p for p in allp if 0 < p < 1]) + 1


def fig_laps():
    lam0 = 3.7
    s0 = float(np.exp(np.interp(lam0, lam, h)))
    N = 6
    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.15),
                            gridspec_kw={"wspace": 0.28})
    x = np.linspace(0, 1, 3000)

    a = axs[0]
    y = g(x, lam0)
    a.plot(x, y, color=BLUE, lw=1.0)
    turns = [0.5]
    edges = np.concatenate([[0.0], turns, [1.0]])
    for i in range(len(edges) - 1):
        if i % 2 == 0:
            a.axvspan(edges[i], edges[i+1], color=RED, alpha=0.07, lw=0)
    a.axvline(0.5, color=RED, lw=0.5, ls="--")
    a.set_title(r"$g_\lambda$:  $\ell_1=2$", fontsize=8)
    a.set_xlabel(r"$x$"); a.set_ylabel(r"$g_\lambda(x)$")

    a = axs[1]
    y = x.copy()
    for _ in range(N):
        y = g(y, lam0)
    a.plot(x, y, color=BLUE, lw=0.6)
    a.set_title(r"$g_\lambda^{%d}$:  $\ell_{%d}=%d$" % (N, N, lap_count(preimages_log(lam0), N)),
                fontsize=8)
    a.set_xlabel(r"$x$")

    a = axs[2]
    y = x.copy()
    for _ in range(N):
        y = np.where(y < 0.5, s0*y, s0*(1 - y))
    a.plot(x, y, color="#2f6b4f", lw=0.6)
    a.set_title(r"$f_s^{%d}$:  $\ell_{%d}=%d$" % (N, N, lap_count(preimages_tent(s0), N)),
                fontsize=8)
    a.set_xlabel(r"$y$")
    for a in axs:
        a.set_xlim(0, 1); a.set_ylim(-0.02, 1.04)
        a.tick_params(labelsize=7); a.grid(alpha=.18, lw=.5)
        a.spines[["top", "right"]].set_visible(False)
    fig.savefig("figs/fig_laps.pdf", bbox_inches="tight")


def fig_model():
    lam0 = 3.7
    s0 = float(np.exp(np.interp(lam0, lam, h)))

    def A_of(x, n=64):
        y = x; w = []
        for _ in range(n):
            w.append(1 if y < 0.5 else -1)
            y = g(y, lam0)
        u = 0.5
        for c in reversed(w):
            u = u/s0 if c > 0 else 1 - u/s0
        return u

    def tent(y):
        return s0*y if y < 0.5 else s0*(1 - y)

    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.15), gridspec_kw={"wspace": 0.28})

    a = axs[0]
    x = np.linspace(0, 1, 400)
    a.plot(x, g(x, lam0), color=BLUE, lw=1.0)
    a.plot([0, 1], [0, 1], color=GREY, lw=0.6, ls=":")
    px, py = [], []
    z = 0.5
    for _ in range(30):
        w = g(z, lam0)
        px += [z, z, w]; py += [z, w, w]; z = w
    a.plot(px, py, color=RED, lw=0.5)
    a.set_title(r"(a)  $g_\lambda$,  $\lambda=3.7$", fontsize=8)
    a.set_xlabel(r"$x$"); a.set_ylabel(r"$g_\lambda(x)$")

    a = axs[1]
    a.plot([0, 0.5, 1], [0, s0/2, 0], color="#2f6b4f", lw=1.0)
    a.plot([0, 1], [0, 1], color=GREY, lw=0.6, ls=":")
    px, py = [], []
    z = 0.5
    for _ in range(30):
        w = tent(z)
        px += [z, z, w]; py += [z, w, w]; z = w
    a.plot(px, py, color=RED, lw=0.5)
    a.set_title(r"(b)  $f_s$,  $s=e^{h(\lambda)}=%.3f$" % s0, fontsize=8)
    a.set_xlabel(r"$y$"); a.set_ylabel(r"$f_s(y)$")

    a = axs[2]
    gx = np.linspace(0, 1, 300)
    a.plot(gx, [A_of(v) for v in gx], color="#333", lw=1.0)
    oL = [0.5]
    for _ in range(6):
        oL.append(g(oL[-1], lam0))
    a.plot(oL, [A_of(v) for v in oL], "o", ms=3.5, color=RED, zorder=5)
    a.plot([0, 1], [0, 1], color=GREY, lw=0.6, ls=":")
    a.set_title(r"(c)  $A$:  $f_s\circ A=A\circ g_\lambda$", fontsize=8)
    a.set_xlabel(r"$x$"); a.set_ylabel(r"$A(x)$")

    for a in axs:
        a.set_xlim(0, 1); a.set_ylim(-0.02, 1.04)
        a.tick_params(labelsize=7); a.grid(alpha=.18, lw=.5)
        a.spines[["top", "right"]].set_visible(False)
    fig.savefig("figs/fig_model.pdf", bbox_inches="tight")


# ---------------------------------------------------------------- figure 4
def eps_seq(l, N):
    x = 0.5
    e = 1.0
    out = np.empty(N + 1)
    out[0] = 1.0
    for k in range(1, N + 1):
        x = g(x, l)
        e = e if x < 0.5 else -e
        out[k] = e
    return out


def root_of(eps, lo=0.45, hi=1 - 1e-12):
    t = np.linspace(lo, hi, 4000)
    D = np.polynomial.polynomial.polyval(t, eps)
    idx = np.where((D[:-1] > 0) & (D[1:] < 0))[0]
    if len(idx) == 0:
        return np.nan
    a, b = t[idx[0]], t[idx[0] + 1]
    for _ in range(80):
        m = 0.5 * (a + b)
        if np.polynomial.polynomial.polyval(m, eps) > 0:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def h_N(l, N):
    t0 = root_of(eps_seq(l, N))
    return -np.log(t0) if np.isfinite(t0) else 0.0


def fig_determinant():
    lam0 = 3.7
    N = 800
    eps = eps_seq(lam0, N)
    tt = np.linspace(0.45, 0.999, 1200)
    D = np.polynomial.polynomial.polyval(tt, eps)
    t0 = root_of(eps)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.35),
                                 gridspec_kw={"width_ratios": [1, 1.15], "wspace": 0.28})
    a1.axvspan(1 - 25 / N, 1.0, color="0.9", lw=0)
    a1.axhline(0, color="0.75", lw=0.6)
    a1.plot(tt, D, color=BLUE, lw=1.0)
    a1.plot([t0], [0], "o", ms=3.5, color=RED, zorder=5)
    a1.annotate(r"$t_0=e^{-h}$", xy=(t0, 0), xytext=(t0 - 0.02, 0.45), fontsize=7,
                color=RED, ha="right", arrowprops=dict(arrowstyle="->", color=RED, lw=0.6))
    a1.text(1 - 25 / N, -0.85, "untrusted", fontsize=6, color=GREY, ha="right", rotation=90,
            va="bottom")
    a1.set_xlim(0.45, 1.0); a1.set_ylim(-1.0, 1.0)
    a1.set_xlabel(r"$t$"); a1.set_ylabel(r"$D_N(t)$")
    a1.set_title(r"(a)  $D_N$ for $\lambda=3.7$, $N=800$", fontsize=8)

    Ns = np.array([16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512])
    for l, col in ((3.575, "#1f77b4"), (3.6, "#ff7f0e"), (3.65, "#2ca02c"), (3.7, "#d62728")):
        ref = h_N(l, 60000)
        xs, es = [], []
        for n in Ns:
            e = abs(h_N(l, int(n)) - ref)
            if e > 3e-15:
                xs.append(ref * n); es.append(e)
        a2.semilogy(xs, es, "o-", ms=2.5, lw=0.8, color=col, label=r"$\lambda=%.3f$" % l)
    a2.axvline(25, color=GREY, lw=0.8, ls="--")
    a2.text(26, 1e-2, r"$hN=25$", fontsize=7, color=GREY)
    a2.set_xlabel(r"$h\,N$"); a2.set_ylabel(r"$|h_N-h|$")
    a2.set_ylim(1e-17, 1e1)
    a2.grid(alpha=0.18, lw=0.5, which="both")
    a2.legend(frameon=False, loc="lower left", handlelength=1.2, fontsize=6)
    a2.set_title(r"(b)  all $\lambda$ collapse onto $e^{-hN}$", fontsize=8)
    fig.savefig("figs/fig_determinant.pdf", bbox_inches="tight")


# ---------------------------------------------------------------- figure 5
def lap_counts(l, nmax):
    S = [0.5]; counts = [1]
    for _ in range(nmax):
        v = []
        for y in S:
            disc = 1 - 4 * y / l
            if disc < 0:
                continue
            r = np.sqrt(disc)
            v.append(0.5 * (1 + r))
            if r > 1e-12:
                v.append(0.5 * (1 - r))
        v = np.sort(np.array(v))
        w = []
        for z in v:
            if not w or z - w[-1] > 1e-11:
                w.append(z)
        S = w; counts.append(len(S))
        if len(S) > 3_000_000:
            break
    return np.array(counts, dtype=float)


def fig_lapgrowth():
    fig, a = plt.subplots(figsize=(6.6, 2.5))
    cols = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    for l, c in zip((3.6, 3.72, 3.83, 3.95, 4.0), cols):
        cnt = lap_counts(l, 24)
        n = np.arange(len(cnt))
        a.plot(n, np.log(cnt), "o", ms=2.2, color=c, label=r"$\lambda=%.2f$" % l)
        hh = h_of(l)
        a.plot(n, np.log(cnt[-1]) + hh * (n - n[-1]), "-", lw=0.7, color=c, alpha=0.7)
    a.set_xlabel(r"$n$")
    a.set_ylabel(r"$\log\,\#\{x: g_\lambda^n(x)=\frac{1}{2}\}$")
    a.grid(alpha=0.18, lw=0.5)
    a.legend(frameon=False, ncol=5, loc="upper left", columnspacing=0.9, handletextpad=0.2)
    a.spines[["top", "right"]].set_visible(False)
    a.set_title(r"straight lines have slope $h(\lambda)$", fontsize=8)
    fig.savefig("figs/fig_lapgrowth.pdf", bbox_inches="tight")


# ---------------------------------------------------------------- figure 6
def density(x0, x1, nc=1600, nr=1000, transient=3000, keep=2000):
    lams = np.linspace(x0, x1, nc)
    x = np.full(nc, 0.5)
    for _ in range(transient):
        x = g(x, lams)
    Z = np.zeros((nr, nc)); idx = np.arange(nc)
    for _ in range(keep):
        x = g(x, lams)
        r = np.clip((x * nr).astype(int), 0, nr - 1)
        np.add.at(Z, (r, idx), 1)
    return np.log1p(Z) / np.log1p(Z.max())


def ink_cmap():
    import matplotlib.colors as mcolors
    return mcolors.LinearSegmentedColormap.from_list(
        "ink", [(1, 1, 1, 0), (0.24, 0.24, 0.24, 0.30), (0.10, 0.10, 0.10, 0.82), (0, 0, 0, 0.95)])


def density_axes(ax, lo, hi, cmap, label=False, top=1.0):
    """The limit-set density of [lo,hi] behind `ax`, on a twin axis."""
    t = ax.twinx()
    t.imshow(density(lo, hi), origin="lower", extent=[lo, hi, 0, 1], aspect="auto",
             cmap=cmap, vmin=0, vmax=1, interpolation="nearest", zorder=1)
    t.set_ylim(0, top)
    if label:
        t.set_yticks([0, 0.5, 1])
        t.set_ylabel(r"$x$   (limit set)", color=GREY, fontsize=7)
    else:
        t.set_yticks([])
    t.tick_params(colors=GREY, length=2, labelsize=6)
    for sp in t.spines.values():
        sp.set_visible(False)
    ax.set_zorder(5)
    ax.patch.set_visible(False)
    return t


# ---------------------------------------------------------------- figure 7
def fig_plateau():
    cmap = ink_cmap()
    panels = [
        (3.820, 3.865, 0.481211825, 3.8284271, 3.8414990, 3.84943, 0.476, 0.497,
         [3.82, 3.84, 3.86], [0.478, 0.482, 0.486, 0.490, 0.494], r"(a)  period 3"),
        (3.902, 3.911, 0.543535072, 3.9055719, 3.9061074, 3.906454, 0.537, 0.551,
         [3.902, 3.906, 3.910], [0.538, 0.542, 0.546, 0.550], r"(b)  period 5"),
        (3.955, 3.967, 0.609377863, 3.9601019, 3.9607687, 3.961198, 0.601, 0.620,
         [3.955, 3.960, 3.965], [0.602, 0.606, 0.610, 0.614, 0.618], r"(c)  period 4"),
    ]
    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.3), gridspec_kw={"wspace": 0.36})
    for k, (a, (lo, hi, h0, sn, pd, acc, ylo, yhi, xt, yt, title)) in enumerate(zip(axs, panels)):
        density_axes(a, lo, hi, cmap, label=(k == 2))
        m = (lam >= lo) & (lam <= hi)
        a.plot(lam[m], h[m], color=BLUE, lw=1.4, zorder=5)
        a.set_xlim(lo, hi); a.set_ylim(ylo, yhi)
        a.set_xticks(xt); a.set_yticks(yt)
        a.set_xlabel(r"$\lambda$", fontsize=8)
        a.set_title(title, fontsize=8)
        a.tick_params(labelsize=6.5)
        a.grid(alpha=0.12, lw=0.5, zorder=0)
        a.axhline(h0, color=GREY, lw=0.6, ls="--")
        for x in (sn, pd):
            a.axvline(x, color="0.4", lw=0.6, ls=":")
        a.axvline(acc, color="0.35", lw=0.5, ls="-.")
        a.spines[["top"]].set_visible(False)
        if k == 0:
            a.set_ylabel(r"$h(\lambda)$", fontsize=8)
    fig.savefig("figs/fig_plateau.pdf", bbox_inches="tight")


if __name__ == "__main__":
    fig_output()
    fig_kneading()
    fig_laps()
    fig_model()
    fig_determinant()
    fig_lapgrowth()
    fig_plateau()
    print("wrote figs/fig_output.pdf, fig_kneading.pdf, fig_laps.pdf, fig_model.pdf,")
    print("      fig_determinant.pdf, fig_lapgrowth.pdf, fig_plateau.pdf")

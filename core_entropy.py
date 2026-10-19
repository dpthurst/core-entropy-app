#!/usr/bin/env python3
"""
core_entropy.py -- core entropy of the logistic family, in Python.

A port of src/CoreEntropy.jl.  For g_λ(x) = λx(1-x) the critical orbit gives
signs s_j = sign g_λ'(x_j) (a return to 1/2 read as the 0-symbol, s_j = -1),
the Milnor-Thurston kneading determinant is D(t) = Σ ε_k t^k with
ε_k = s_1···s_k, and the core entropy is h(λ) = log(1/t₀), where t₀ is the
smallest zero of D in (0,1); h = 0 if there is none.

Truncating at k = N leaves the tail bound |D - D_N| ≤ t^(N+1)/(1-t), so a zero
is accepted only once h·N ≥ 25; otherwise N is doubled, up to 51200.  The sweep
is vectorised over λ: one occupancy-free Horner pass per t, all parameters at
once.

    python3 core_entropy.py          # writes core_entropy_data.csv
"""
import numpy as np

LOG2 = np.log(2.0)
LAM_INF = 3.5699456718        # Feigenbaum accumulation
LAM_P3 = 3.8318740553         # period-3 centre
RESOLUTION = 25.0


def lambda_grid():
    """The non-uniform grid: dense where h moves, sparse on the flat part."""
    pts = np.concatenate([
        np.linspace(2.00, 3.40, 600),
        np.linspace(3.40, 3.55, 600),
        np.linspace(3.55, 3.58, 1500),
        np.linspace(3.58, 3.72, 4000),
        np.linspace(3.72, 4.00, 22000)])
    return np.unique(pts)


def sign_matrix(lams, N):
    """eps[i, k] = ε_k for lams[i], k = 0..N (int8, entries ±1)."""
    n = len(lams)
    eps = np.empty((n, N + 1), dtype=np.int8)
    eps[:, 0] = 1
    e = np.ones(n, dtype=np.int8)
    x = np.full(n, 0.5)
    for k in range(1, N + 1):
        x = lams * x * (1.0 - x)
        e = np.where(x < 0.5, e, -e).astype(np.int8)
        eps[:, k] = e
    return eps


def _D_rows(eps, t):
    """D_N(t_i) for each row i, by Horner in t over the coefficient rows."""
    s = np.zeros(len(t))
    for k in range(eps.shape[1] - 1, -1, -1):
        s = s * t + eps[:, k]
    return s


def _D_grid(eps, ts):
    """D_N at every (parameter, t) pair."""
    out = np.empty((eps.shape[0], len(ts)))
    for j, t in enumerate(ts):
        out[:, j] = _D_rows(eps, np.full(eps.shape[0], t))
    return out


def core_entropy_one(lam, N0=800, Nmax=51200, resolution=RESOLUTION):
    """Scalar version with the adaptive order; used for the small-h parameters."""
    N = N0
    tg = None
    while True:
        eps = sign_matrix(np.array([lam]), N)[0]
        tmax = 1.0 - resolution / N
        t0 = _smallest_root(eps, tg, tmax)
        if t0 is None:
            if N < Nmax and resolution / N > 0.008:
                N = min(4 * N, Nmax)
                continue
            return 0.0
        h = -np.log(t0)
        if h * N >= resolution or N >= Nmax:
            return h
        tg = t0
        N = min(2 * N, Nmax)


def _smallest_root(eps, tguess, tmax):
    """First sign change of D from + to - on (0.45, tmax), by bisection."""
    def D(t):
        s = 0.0
        for k in range(len(eps) - 1, -1, -1):
            s = s * t + eps[k]
        return s
    ts = np.linspace(0.45, tmax, 80)
    prev, fprev = 0.45, D(0.45)
    for t in ts[1:]:
        ft = D(t)
        if fprev > 0 > ft:
            a, b = prev, t
            for _ in range(70):
                m = 0.5 * (a + b)
                if D(m) > 0:
                    a = m
                else:
                    b = m
            return 0.5 * (a + b)
        prev, fprev = t, ft
    return None


def compute_curve(lams, N=800, verbose=True):
    """h(λ) for every λ on the grid: vectorised sweep, then scalar refinement."""
    n = len(lams)
    eps = sign_matrix(lams, N)
    tmax = 1.0 - RESOLUTION / N
    ts = np.linspace(0.45, tmax, 96)
    D = _D_grid(eps, ts)

    cross = (D[:, :-1] > 0) & (D[:, 1:] < 0)
    found = cross.any(axis=1)
    j = np.argmax(cross, axis=1)
    lo = np.where(found, ts[np.clip(j, 0, len(ts) - 2)], 0.45)
    hi = np.where(found, ts[np.clip(j + 1, 1, len(ts) - 1)], 0.45)
    for _ in range(48):
        mid = 0.5 * (lo + hi)
        up = _D_rows(eps, mid) > 0
        lo = np.where(up, mid, lo)
        hi = np.where(up, hi, mid)
    t0 = 0.5 * (lo + hi)

    h = np.where(found, -np.log(np.where(found, t0, 0.5)), 0.0)

    # the vectorised pass used N = 800; anything it could not resolve gets the
    # scalar treatment (h·N < 25, i.e. h < 0.031), and so does every parameter
    # above λ∞ that came out as zero
    todo = [(h > 0) & (h * N < RESOLUTION)] if verbose else []
    if verbose:
        idx = np.where((h > 0) & (h * N < RESOLUTION))[0]
        print(f"   {n} parameters, {np.count_nonzero(found)} zeros found, "
              f"{len(idx)} need refining")
        for i in idx:
            h[i] = core_entropy_one(lams[i])
        idx0 = np.where((h == 0) & (lams > 3.55))[0]
        print(f"   {len(idx0)} zero-entropy candidates above 3.55, rechecking")
        for i in idx0:
            h[i] = core_entropy_one(lams[i])
    return h


def h_at(lams, hs, lam):
    """Value of h at the grid parameter nearest to λ."""
    i = int(np.clip(np.searchsorted(lams, lam), 0, len(lams) - 1))
    if i > 0 and abs(lams[i - 1] - lam) < abs(lams[i] - lam):
        i -= 1
    return hs[i]


def write_csv(path, lams, hs):
    c = lams / 2 - lams ** 2 / 4
    np.savetxt(path, np.column_stack([lams, c, hs]), delimiter=",", fmt="%.12f",
               header="lambda,c,h", comments="")
    return path


if __name__ == "__main__":
    import os, time
    t = time.time()
    lams = lambda_grid()
    print(f"computing h(λ) on {len(lams)} parameters ...")
    hs = compute_curve(lams)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "core_entropy_data.csv")
    write_csv(out, lams, hs)
    print(f"wrote {out} in {time.time()-t:.1f}s")
    print(f"h(2) = {hs[0]:.1e}   h(4) = {hs[-1]:.12f}   (log 2 = {np.log(2):.12f})")

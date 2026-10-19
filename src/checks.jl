# ================================================================ checks

"""
    preimage_growth(λ, nmax) -> Vector{Int}

Number of x ∈ [0,1] with g_λ^n(x) = 1/2 for n = 0..nmax -- the number of laps
of g_λ^n, which grows like e^{hn}.  Uses the dynamics only (no kneading theory),
so it is an independent check of h.
"""
function preimage_growth(λ::Float64, nmax::Int; cap::Int = 4_000_000)
    S = [0.5]
    counts = Int[1]
    for _ in 1:nmax
        v = Float64[]
        for y in S
            disc = 1.0 - 4.0 * y / λ
            disc < 0 && continue
            r = sqrt(disc)
            push!(v, 0.5 * (1 + r))
            r > 1e-12 && push!(v, 0.5 * (1 - r))
        end
        sort!(v)
        w = Float64[]
        for z in v                          # merge numerically coincident preimages
            (isempty(w) || z - w[end] > 1e-11) && push!(w, z)
        end
        S = w
        push!(counts, length(S))
        length(S) > cap && break
    end
    return counts
end

function run_checks()
    ok = true
    println("core entropy: validation")
    println("=======================")
    println()
    println("1. parameters whose entropy is known exactly")
    c3 = -1.7548776662466927                 # centre of the period-3 hyperbolic component
    lam3 = 1 + sqrt(1 - 4c3)
    tests = [
        (2.0,           0.0,      "critical point fixed (c = 0)"),
        (3.0,           0.0,      "parabolic fixed point (c = -3/4)"),
        (1 + sqrt(5.0), 0.0,      "basilica centre (c = -1)"),
        (lam3,          log(PHI), "period-3 centre (c = $(round(c3, digits = 7))): log τ"),
        (4.0,           LOG2,     "Chebyshev map (c = -2): log 2"),
    ]
    for (λ, want, name) in tests
        h, _, N = core_entropy(λ; N0 = 1600, Nmax = 102400, resolution = 40.0)
        err = abs(h - want)
        ok &= err < 1e-9
        @printf("   h(%.10f) = %.12f   expected %.12f   |err| = %.1e   N = %d\n      %s\n",
                λ, h, want, err, N, name)
    end

    println()
    println("2. independent check: growth rate of the number of preimages of the")
    println("   critical point (laps of g^n); finite-n estimates converge slowly")
    println("   when h is small, hence the 1% tolerance")
    for λ in (3.62, 3.72, 3.83, 3.95, 4.0)
        h, _, _ = core_entropy(λ; N0 = 1600, resolution = 40.0)
        cnt = preimage_growth(λ, 26)
        n = length(cnt) - 1
        hnum = (log(cnt[end]) - log(cnt[max(1, end - 4)])) / 4
        err = abs(hnum - h) / max(h, 1e-9)
        ok &= err < 0.01
        @printf("   λ = %.4f   h = %.9f   laps(n=%d) = %-9d  growth estimate %.6f   rel. diff %.1e\n",
                λ, h, n, cnt[end], hnum, err)
    end

    println()
    println("3. monotonicity on a sweep of the whole interval")
    lams = collect(range(2.0, 4.0; length = 4001))
    hs = [core_entropy(l)[1] for l in lams]
    dmin = minimum(diff(hs))
    ok &= dmin > -1e-9
    @printf("   min h(λ+Δ) - h(λ) over 4001 points = %.2e  (must be ≥ 0)\n", dmin)
    @printf("   h(2.0) = %.1e,  h(4.0) = %.12f,  log 2 = %.12f\n", hs[1], hs[end], LOG2)

    println()
    println(ok ? "all checks passed." : "SOME CHECKS FAILED.")
    return ok
end


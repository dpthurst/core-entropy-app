# ================================================================ checks

"""
    preimage_growth(λ, nmax) -> Vector{Int}

Number of x ∈ [0,1] with g_λ^n(x) = 1/2 for n = 0..nmax.  It grows at the
same exponential rate e^{hn} as the number of laps of g_λ^n, and uses the
dynamics only (no kneading theory), so it is an independent check of h.
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
        # Merge numerically coincident preimages.  The tolerance has to stay far below
        # the smallest gap, which falls like λ^-n: at 1e-11 distinct preimages were being
        # merged from n = 21, and λ = 3.62 counted 3068 instead of 3106 at n = 26.  At
        # 1e-15 the counts agree with a 35-digit computation as far as double precision
        # reaches.
        for z in v
            (isempty(w) || z - w[end] > 1e-15) && push!(w, z)
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
        (lam3,          log(TAU), "period-3 centre (c = $(round(c3, digits = 7))): log τ"),
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
    println("   when h is small, hence the 2% tolerance")
    for λ in (3.62, 3.72, 3.83, 3.95, 4.0)
        h, _, _ = core_entropy(λ; N0 = 1600, resolution = 40.0)
        cnt = preimage_growth(λ, 26)
        n = length(cnt) - 1
        hnum = (log(cnt[end]) - log(cnt[max(1, end - 4)])) / 4
        err = abs(hnum - h) / max(h, 1e-9)
        # The finite-n estimate is biased high by the subexponential correction, and the
        # bias is worst where h is small: at lambda = 3.62 it is 1.1% at n = 26, which is
        # why the threshold is 2% rather than 1%.  The 1% of earlier versions was only met
        # because merged preimages made the counts short.
        ok &= err < 0.02
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
    println("4. the critical itinerary is an admissible kneading sequence")
    println("   a kneading sequence must dominate every shift of itself in the kneading")
    println("   order: at the first difference, 0 < 1 when the number of 1s before it is")
    println("   even and 1 < 0 when it is odd; parameters where the orbit returns are skipped")
    tested = 0
    violations = 0
    firstbad = (NaN, 0)
    for λ in range(3.58, 4.0; length = 400)
        s = critical_itinerary(λ)
        occursin('*', s) && continue
        tested += 1
        for k in 2:length(s)
            if !kneading_ge(s, s[k:end])
                violations += 1
                isnan(firstbad[1]) && (firstbad = (λ, k - 1))
                break
            end
        end
    end
    ok &= violations == 0
    @printf("   %d parameters tested, %d violations", tested, violations)
    isnan(firstbad[1]) || @printf("   (first at λ = %.6f, shift %d)", firstbad[1], firstbad[2])
    println()

    println()
    println(ok ? "all checks passed." : "SOME CHECKS FAILED.")
    return ok
end


"""
The itinerary of the critical point: '0' left of 1/2, '1' right, '*' on a return.
"""
function critical_itinerary(λ::Float64, n::Int = 48; tol::Float64 = 1e-6)
    out = Vector{Char}(undef, 0)
    x = 0.5
    for _ in 1:n
        x = λ * x * (1 - x)
        push!(out, abs(x - 0.5) < tol ? '*' : (x < 0.5 ? '0' : '1'))
    end
    return String(out)
end

"""
    kneading_ge(u, v) -> Bool

Whether u ≥ v in the kneading order: at the first place they differ, 0 < 1 when the
number of 1s in u before that place is even, and 1 < 0 when it is odd.
"""
function kneading_ge(u::AbstractString, v::AbstractString)
    for i in 1:min(length(u), length(v))
        a, b = u[i], v[i]
        a == b && continue
        even = iseven(count(==('1'), u[1:i-1]))
        return a != (even ? '0' : '1')
    end
    return true
end


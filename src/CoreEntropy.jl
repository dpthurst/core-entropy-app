"""
    CoreEntropy

The computational core of the core-entropy plot, shared by the two front ends
(`core_entropy.jl`, which writes the interactive HTML page, and
`core_entropy_makie.jl`, which opens a Makie window).  Nothing here depends on a
plotting library.

Core entropy of the real quadratic family g_λ(x) = λx(1-x), λ ∈ [2,4].
g_λ is affinely conjugate to z ↦ z² + c with c = λ/2 - λ²/4 ∈ [-2,0]; for real c
the Hubbard tree is a real interval and the core entropy is the topological
entropy of g_λ on [0,1].  It is a continuous, nondecreasing devil's staircase from
h(2) = 0 to h(4) = log 2 (Tiozzo, Adv. Math. 273 (2015)).

Computation.  With x_j = g_λ^j(1/2) and s_j = sign g_λ'(x_j) (a return
x_j = 1/2 read as the 0-symbol, s_j = -1), the Milnor-Thurston kneading
determinant is D(t) = Σ ε_k t^k with ε_k = ∏_{j≤k} s_j, and h = log(1/t₀) where
t₀ is the smallest zero of D in (0,1): 1/t₀ is the growth rate of the number of
laps of g^n.  No zero means h = 0.  Truncating at k = N leaves the tail bound
|D - D_N| ≤ t^(N+1)/(1-t), so a zero is accepted only once (1-t₀)N ≥ 25.
"""
module CoreEntropy

using Printf

# landmarks of the family (λ ↦ c = λ/2 − λ²/4)
const LOG2       = log(2.0)
const PHI        = (1 + sqrt(5.0)) / 2
const LAMBDA_INF = 3.5699456718     # Feigenbaum accumulation, c ≈ -1.4011552
const LAMBDA_P3  = 3.8318740553     # period-3 centre,       c ≈ -1.7548777

export LOG2, PHI, LAMBDA_INF, LAMBDA_P3
export kneading_eps, kneading_D, smallest_root, core_entropy, preimage_growth
export lambda_grid, compute_curve, write_csv, h_at, critical_orbit, run_checks

# ================================================================ kneading theory

"""
    kneading_eps(λ, N) -> Vector{Float64}

ε₀,…,ε_N with ε_k = ∏_{j=1..k} s_j and s_j = sign g_λ'(g_λ^j(1/2)).
"""
function kneading_eps(λ::Float64, N::Int)
    eps = Vector{Float64}(undef, N + 1)
    eps[1] = 1.0
    e = 1.0
    x = 0.5
    @inbounds for k in 1:N
        x = λ * x * (1.0 - x)
        e = x < 0.5 ? e : -e
        eps[k+1] = e
    end
    return eps
end

"""Kneading determinant D(t) truncated at N, evaluated by Horner's rule."""
@inline function kneading_D(eps::Vector{Float64}, t::Float64)
    s = 0.0
    @inbounds for k in length(eps):-1:1
        s = muladd(s, t, eps[k])
    end
    return s
end

function _bisect(eps::Vector{Float64}, a::Float64, b::Float64)
    for _ in 1:80
        m = 0.5 * (a + b)
        (m <= a || m >= b) && break
        if kneading_D(eps, m) > 0
            a = m
        else
            b = m
        end
        (b - a) <= 4e-16 * max(a, b) && break
    end
    return 0.5 * (a + b)
end

"""
    smallest_root(eps; tguess = NaN, tmax = 1 - 1e-12)

Smallest t in (0.45, tmax) at which the truncated determinant changes sign from
+ to -.  Returns NaN if there is none.  `tguess` (the zero found at a nearby
parameter) is tried first as a bracket, which makes a sweep over a grid fast.
"""
function smallest_root(eps::Vector{Float64}; tguess::Float64 = NaN,
                       tmax::Float64 = 1 - 1e-12)
    if !isnan(tguess)                       # bracket around the previous zero
        δ = 1e-4
        for _ in 1:14
            a = max(tguess - δ, 0.45)
            b = min(tguess + δ, tmax)
            if a < b && kneading_D(eps, a) > 0 && kneading_D(eps, b) < 0
                return _bisect(eps, a, b)
            end
            δ *= 4
        end
    end
    tprev, fprev = 0.45, kneading_D(eps, 0.45)   # D(0)=1>0, D stays >0 near 0
    M = 80
    for j in 1:M
        t = 0.45 + (tmax - 0.45) * j / M
        ft = kneading_D(eps, t)
        if fprev > 0 && ft < 0
            return _bisect(eps, tprev, t)
        end
        tprev, fprev = t, ft
    end
    return NaN
end

"""
    core_entropy(λ; N0=800, Nmax=51200, tguess=NaN, resolution=25.0) -> (h, t₀, N)

Core entropy h(λ) = log(1/t₀), its zero t₀, and the number of kneading terms
used.  A zero is trusted only once (1-t₀)N ≥ resolution, otherwise N is doubled;
if no zero is found the entropy is reported as 0, resolved down to h ≈
resolution/N.
"""
function core_entropy(λ::Float64; N0::Int = 800, Nmax::Int = 51200,
                      tguess::Float64 = NaN, resolution::Float64 = 25.0)
    N = N0
    tg = tguess
    while true
        eps = kneading_eps(λ, N)
        tmax = 1.0 - resolution / N
        t0 = smallest_root(eps; tguess = tg, tmax = tmax)
        if isnan(t0)                       # nothing resolvable at this N:
            if N < Nmax && resolution / N > 0.008   # one more look before saying 0
                N = min(4N, Nmax)
                continue
            end
            return (0.0, NaN, N)
        end
        h = -log(t0)
        if h * N >= resolution || N >= Nmax
            return (h, t0, N)
        end
        tg = t0
        N = min(2N, Nmax)
    end
end

# ================================================================ data

"""
    lambda_grid() -> Vector{Float64}

Non-uniform grid: dense where h moves (λ ≳ 3.57) and sparse on the flat part.
"""
function lambda_grid()
    pts = Float64[]
    append!(pts, range(2.00, 3.40; length = 600))
    append!(pts, range(3.40, 3.55; length = 600))
    append!(pts, range(3.55, 3.58; length = 1500))
    append!(pts, range(3.58, 3.72; length = 4000))
    append!(pts, range(3.72, 4.00; length = 22000))
    sort!(unique!(pts))
    return pts
end

function compute_curve(lams::Vector{Float64}; verbose::Bool = true)
    n = length(lams)
    hs = zeros(n)
    t0s = fill(NaN, n)
    ns = zeros(Int, n)
    tg = NaN
    tstart = time()
    for i in 1:n
        h, t0, N = core_entropy(lams[i]; tguess = tg)
        hs[i] = h
        t0s[i] = t0
        ns[i] = N
        tg = isnan(t0) ? NaN : t0
        if verbose && i % 5000 == 0
            @printf("   %5d / %d   (%.1f s)\n", i, n, time() - tstart)
        end
    end
    verbose && @printf("   %d parameters in %.1f s, largest N used: %d\n",
                       n, time() - tstart, maximum(ns))
    return hs, t0s, ns
end

function write_csv(path::String, lams, hs)
    open(path, "w") do io
        println(io, "lambda,c,h")
        for (l, h) in zip(lams, hs)
            @printf(io, "%.12f,%.12f,%.12f\n", l, l / 2 - l^2 / 4, h)
        end
    end
    return path
end


# ---------------------------------------------------------------- shared utilities

"""
    h_at(lams, hs, λ) -> Float64

Value of `h` at the grid parameter nearest to `λ` (the grid is not uniform).
"""
function h_at(lams::AbstractVector{Float64}, hs::AbstractVector{Float64}, λ::Real)
    i = clamp(searchsortedlast(lams, λ), 1, length(lams) - 1)
    return λ - lams[i] <= lams[i+1] - λ ? hs[i] : hs[i+1]
end

"""
    critical_orbit(λ; n = 140) -> Vector{NTuple{2,Float64}}

The first `n` points `(x_k, g_λ(x_k))` of the critical orbit, starting at `x_0 = 1/2`.
"""
function critical_orbit(λ::Real; n::Int = 140)
    pts = NTuple{2,Float64}[]
    x = 0.5
    for _ in 1:n
        y = λ * x * (1 - x)
        push!(pts, (x, y))
        x = y
    end
    return pts
end

include("checks.jl")

end # module CoreEntropy

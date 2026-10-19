#!/usr/bin/env julia
#
# core_entropy.jl -- core entropy of the real quadratic family, HTML front end
#
#   g_λ(x) = λ x (1 - x),     2 ≤ λ ≤ 4,   x ∈ [0,1]
#
# Computes h(λ), the topological entropy of g_λ on its core (see
# src/CoreEntropy.jl for the mathematics and the algorithm), and writes
#
#   core_entropy_data.csv    lambda,c,h
#   core_entropy.html        self-contained interactive plot
#
# Usage
#   julia core_entropy.jl              # data + interactive HTML
#   julia core_entropy.jl --check      # validation against exactly known values
#   julia core_entropy.jl --cdn        # smaller HTML, plotly.js from the CDN
#   julia core_entropy.jl --help
#
# Layout
#   src/CoreEntropy.jl     kneading theory, entropy solver, grid, CSV, helpers
#   src/checks.jl          the validation suite (CoreEntropy.run_checks)
#   src/HTMLPage.jl        the Plotly page
#   core_entropy_makie.jl  a second front end over the same core

using Printf

include(joinpath(@__DIR__, "src", "CoreEntropy.jl"))
include(joinpath(@__DIR__, "src", "HTMLPage.jl"))

using .CoreEntropy
using .HTMLPage: write_html

function main(args)
    if "--help" in args || "-h" in args
        println(read(joinpath(@__DIR__, "README.md"), String))
        return
    end
    if "--check" in args
        ok = run_checks()
        exit(ok ? 0 : 1)
    end

    lams = lambda_grid()
    println("core entropy of g_λ(x) = λx(1-x),  λ ∈ [2,4]")
    println("computing $(length(lams)) parameters ...")
    hs, _, _ = compute_curve(lams)

    csv = write_csv(joinpath(@__DIR__, "core_entropy_data.csv"), lams, hs)
    println("wrote $(csv)")
    html = write_html(joinpath(@__DIR__, "core_entropy.html"), lams, hs;
                      cdn = "--cdn" in args)
    @printf("wrote %s (%.1f MB)\n", html, filesize(html) / 2^20)
    @printf("h(2) = %.1e,  h(4) = %.12f  (log 2 = %.12f)\n", hs[1], hs[end], LOG2)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

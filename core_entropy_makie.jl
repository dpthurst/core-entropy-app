#!/usr/bin/env julia
#
# core_entropy_makie.jl -- native Makie window for the core entropy h(λ)
#
# DEPRECATED.  Kept for reference; not maintained.  Use tent_map_app.py, which is
# the same window in Python and matplotlib and needs no Julia, no Makie and no
# display to render a frame:
#
#     python3 tent_map_app.py                  # the window
#     python3 tent_map_app.py --save f.png --exact
#
# Until recently it did not even parse: the Axis(...) call below was missing its
# closing parenthesis, which Julia reports as "Expected `)` or `,`" at the end of
# the xlabel/ylabel line -- on 1.11 and 1.12 alike, so it was never a version
# difference.  That is fixed and the file parses cleanly under Julia 1.11.3 and
# 1.12.5.  What remains untested is the drawing, which needs Makie and a display;
# that is why this stays deprecated in favour of tent_map_app.py.
#
# ---- the original description, for reference -------------------------------
#
# A second front end over src/CoreEntropy.jl, which does all of the computation;
# this file only draws.  Requires Makie and GLMakie:
#     julia -e 'using Pkg; Pkg.add(["Makie", "GLMakie"])'
#
#     julia core_entropy_makie.jl
#
# Drag the slider to move λ: the marker follows the curve and the right panel
# redraws g_λ on [0,1] with its critical orbit.  Pan and zoom as usual.
#
# NOT TESTED: written in an environment with no display and no Makie, so the
# plotting calls below have never been executed.  The core it calls is tested
# (`julia core_entropy.jl --check`).

using GLMakie

include(joinpath(@__DIR__, "src", "CoreEntropy.jl"))
using .CoreEntropy

function makie_view()
    λ0 = LAMBDA_P3
    lams = lambda_grid()
    hs, _, _ = compute_curve(lams; verbose = false)

    fig = Figure(size = (1180, 560), fontsize = 13)

    ax1 = Axis(fig[1:2, 1];
               title = "core entropy of  g_λ(x) = λx(1−x)",
               xlabel = "λ", ylabel = "h(λ)")
    ax2 = Axis(fig[1, 2]; aspect = DataAspect(),
               title = "g_λ on [0,1]", xlabel = "x", ylabel = "g_λ(x)",
               limits = ((0.0, 1.0), (0.0, 1.02)))
    colsize!(fig.layout, 1, Relative(0.62))

    sl = Slider(fig[2, 2]; range = 2.0:0.001:4.0, startvalue = λ0, tellwidth = true)

    # ---- the entropy curve
    lines!(ax1, lams, hs; color = "#1c5f9e", linewidth = 1.5)
    vlines!(ax1, [LAMBDA_INF, LAMBDA_P3]; color = (:gray, 0.5), linestyle = :dash)
    hlines!(ax1, [LOG2]; color = (:gray, 0.3), linestyle = :dot)
    xlims!(ax1, 2.0, 4.0)
    ylims!(ax1, -0.02, 0.70)
    text!(ax1, 3.575, 0.655; text = "λ∞", fontsize = 11, color = (:gray, 0.8))
    text!(ax1, 3.79, 0.655; text = "period 3", fontsize = 11, color = (:gray, 0.8))

    marker = @lift([Point2f($sl.value, h_at(lams, hs, $sl.value))])
    scatter!(ax1, marker; color = "#c1272d", markersize = 13,
             strokewidth = 1, strokecolor = :white)

    # ---- the map
    xs = range(0.0, 1.0; length = 400)
    ys = Observable(λ0 .* xs .* (1 .- xs))
    lines!(ax2, xs, ys; color = "#1c5f9e", linewidth = 2)
    lines!(ax2, [0.0, 1.0], [0.0, 1.0]; color = (:gray, 0.6), linestyle = :dot)
    vlines!(ax2, [0.5]; color = (:gray, 0.4), linestyle = :dash)
    as_points(λ) = [Point2f(p[1], p[2]) for p in critical_orbit(λ)]
    orb = Observable(as_points(λ0))
    scatter!(ax2, orb; color = ("#e08a1e", 0.75), markersize = 6)

    on(sl.value) do λ
        ys[] = λ .* xs .* (1 .- xs)
        orb[] = as_points(λ)
    end

    # ---- the readout
    info = lift(sl.value) do λ
        h = h_at(lams, hs, λ)
        c = λ / 2 - λ^2 / 4
        string("λ = ", round(λ, digits = 6), "      c = ", round(c, digits = 6),
               "\nh = ", round(h, digits = 6), "      h / log 2 = ",
               round(h / LOG2, digits = 4))
    end
    Label(fig[3, 2], info; halign = :left, tellwidth = false, fontsize = 12)

    resize_to_layout!(fig)
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    fig = makie_view()
    screen = display(GLMakie.Screen(), fig)
    wait(screen)
end

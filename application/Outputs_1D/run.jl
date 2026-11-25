# application/Outputs_1D/run.jl -- derived quantities ("cooking recipes") with OPT precision, 1-D test
#
# From a terminal, in the repository root (on stella: inside srun / sbatch, never on stella*):
#   julia --project=. application/Outputs_1D/run.jl                 # one grid, G = 8 points per wavelength
#   julia --project=. application/Outputs_1D/run.jl compare G=5
#   julia --project=. application/Outputs_1D/run.jl convergence     # error vs points per wavelength, all methods
#   julia --project=. application/Outputs_1D/run.jl convergence schemes=OPT3,OPT5
# From a REPL:
#   julia> include("application/Outputs_1D/run.jl")
#   julia> r = outputs1DCompare(G=6); r.err
#
# Outputs: application/Outputs_1D/outputs/<task>/ (PNG + CSV)

FLEXOPT_MODULES = [:flexOPT]
include(joinpath(@__DIR__, "..", "common", "appBootstrap.jl"))
using CairoMakie, DelimitedFiles
CairoMakie.activate!()
include(joinpath(@__DIR__, "app.jl"))

outputs1D(args...) = (global APP_ARGS = collect(String, string.(args)); include(joinpath(@__DIR__, "run.jl")))

OUTPUTS1D_COLORS = Dict("FD3 (Fornberg)" => :gray60, "FD5 (Fornberg)" => :black,
    "convFD3 lumped" => :orange, "convFD3 consistent" => :darkorange3,
    "OPT3 lumped" => :dodgerblue, "OPT3 consistent" => :navy,
    "OPT5 lumped" => :tomato, "OPT5 consistent" => :darkred)
_o1color(m) = get(OUTPUTS1D_COLORS, m, :purple)

let A = appArgs(@isdefined(APP_ARGS) ? APP_ARGS : ARGS; G=8.0, L=60.0, schemes=["convFD3", "OPT3", "OPT5"])
    task = isempty(A.positional) ? "compare" : A.positional[1]
    schemes = A.schemes isa AbstractVector ? String.(string.(A.schemes)) : [string(A.schemes)]
    println("Outputs_1D: task=$task  ", A)

    if task == "compare"
        r = outputs1DCompare(; G=Float64(A.G), L=Float64(A.L), schemes=schemes)
        dir = appOutputDir(@__DIR__, "compare")
        methods = sort(collect(keys(r.out)))
        fig = Figure(size=(1100, 800))
        ax1 = Axis(fig[1, 1]; ylabel="q = a ∂x u",
                   title=@sprintf("q = a(x) ∂x u,  u = cos x,  a = 2 + sin(x/3),  G = %.1f points per wavelength", r.G))
        lines!(ax1, r.X, r.q; color=:black, linewidth=3, label="exact")
        for m in methods
            scatterlines!(ax1, r.X, r.out[m]; color=_o1color(m), markersize=4, linewidth=1, label=m)
        end
        axislegend(ax1; position=:rb, nbanks=3, labelsize=10)
        ax2 = Axis(fig[2, 1]; xlabel="x", ylabel="|error| (log)", yscale=log10)
        for m in methods
            e = abs.(r.out[m] .- r.q) .+ 1e-16
            lines!(ax2, r.X, e; color=_o1color(m), label=@sprintf("%s  (%.1e)", m, r.err[m]))
        end
        vspan!(ax2, [r.X[1], r.X[last(r.interior)]], [r.X[first(r.interior)], r.X[end]]; color=(:gray, 0.15))
        axislegend(ax2; position=:rb, nbanks=2, labelsize=10)
        appSave(fig, dir, @sprintf("compare_G%.1f", r.G))
        writedlm(joinpath(dir, @sprintf("compare_G%.1f.csv", r.G)),
                 vcat(permutedims(vcat(["x", "exact"], methods)), hcat(r.X, r.q, [r.out[m] for m in methods]...)), ',')
        global outputs1DLast = r

    elseif task == "convergence"
        c = outputs1DConvergence(; L=Float64(A.L), schemes=schemes)
        dir = appOutputDir(@__DIR__, "convergence")
        fig = Figure(size=(900, 650))
        ax = Axis(fig[1, 1]; xlabel="points per wavelength G", ylabel="relative L2 error of q = a ∂x u",
                  xscale=log10, yscale=log10, title="Derived quantity: conventional FD vs OPT cooking recipes")
        for m in c.methods
            scatterlines!(ax, c.Gs, c.err[m]; color=_o1color(m), label=m,
                          linestyle=occursin("lumped", m) ? :dash : :solid)
        end
        for (p, y0) in ((2, 0.5), (4, 0.5))
            lines!(ax, c.Gs, y0 .* (c.Gs[1] ./ c.Gs) .^ p; color=:gray80, linestyle=:dot)
        end
        axislegend(ax; position=:lb, labelsize=11)
        appSave(fig, dir, "convergence")
        writedlm(joinpath(dir, "convergence.csv"),
                 vcat(permutedims(vcat(["G"], c.methods)), hcat(c.Gs, [c.err[m] for m in c.methods]...)), ',')
        global outputs1DConvergenceLast = c
    else
        error("unknown task '$task' (compare, convergence)")
    end
end

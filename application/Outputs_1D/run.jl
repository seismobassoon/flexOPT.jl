# application/Outputs_1D/run.jl -- derived quantities chosen as PDE expressions, cooked with OPT recipes (1-D test)
#
# From a terminal, in the repository root (on stella: inside srun / sbatch, never on stella*):
#   julia --project=. application/Outputs_1D/run.jl                          # compare, output "flux" (a ∂x u), G = 8
#   julia --project=. application/Outputs_1D/run.jl compare name=mixed G=6   # a ∂x u + b u
#   julia --project=. application/Outputs_1D/run.jl convergence              # error vs points per wavelength
#   julia --project=. application/Outputs_1D/run.jl custom expr="a*∂x(∂x(u))+b*∂x(u)"   # your own expression
#        (materials a, b, c = 2+sin(x/3), 1+cos(x/5)/2, 1.5; field u = cos x; no commas in expr)
# From a REPL (one line per statement):
#   julia> include("application/Outputs_1D/run.jl")
#   julia> registerOutput!("stress", (u, μ) -> μ * ∂x(u))
#   julia> r = outputs1DCompare("stress"; G=6); r.err
#   julia> lumped, consistent = cookOutput("stress", uvals, [μvals], Δx; scheme="OPT3", mass="Gal3")
#
# Outputs: application/Outputs_1D/outputs/<task>/ (PNG + CSV)

FLEXOPT_MODULES = [:flexOPT]
include(joinpath(@__DIR__, "..", "common", "appBootstrap.jl"))
using CairoMakie, DelimitedFiles
CairoMakie.activate!()
include(joinpath(@__DIR__, "app.jl"))

outputs1D(args...) = (global APP_ARGS = collect(String, string.(args)); include(joinpath(@__DIR__, "run.jl")))

OUTPUTS1D_ARGS = appArgs(@isdefined(APP_ARGS) ? APP_ARGS : ARGS; name="flux", G=8.0, L=60.0, expr="")
OUTPUTS1D_TASK = isempty(OUTPUTS1D_ARGS.positional) ? "compare" : OUTPUTS1D_ARGS.positional[1]
println("Outputs_1D: task=$OUTPUTS1D_TASK  ", OUTPUTS1D_ARGS)

# a user expression from the command line: registered here, at top level (world age), as output "custom"
if OUTPUTS1D_TASK == "custom"
    _e = string(OUTPUTS1D_ARGS.expr)
    isempty(_e) && error("custom needs expr=\"...\" (field u, materials a, b, c)")
    _m = [m for m in ("a", "b", "c") if occursin(Regex("\\b$m\\b"), _e)]
    registerOutput!("custom", eval(Meta.parse("(u" * join(", " .* _m) * ") -> " * _e)))
end

_o1color(m) = occursin("Gal", m) ? (occursin("OPT5", m) ? :darkred : occursin("OPT3", m) ? :navy : :seagreen) :
              occursin("OPT5", m) ? :tomato : occursin("OPT3", m) ? :dodgerblue :
              occursin("convFD5", m) ? :black : :gray55

let A = OUTPUTS1D_ARGS, task = OUTPUTS1D_TASK
    name = task == "custom" ? "custom" : string(A.name)
    if task in ("compare", "custom")
        r = outputs1DCompare(name; G=Float64(A.G), L=Float64(A.L))
        dir = appOutputDir(@__DIR__, task)
        methods = sort(collect(keys(r.out)); by=k -> r.err[k])
        fig = Figure(size=(1150, 850))
        ax1 = Axis(fig[1, 1]; ylabel="q", title=@sprintf("q = %s,  u = cos x,  G = %.1f points per wavelength", string(OUTPUT_REGISTRY[name].expr), r.G))
        lines!(ax1, r.X, r.q; color=:black, linewidth=3, label="exact")
        for m in methods
            occursin("consistent", m) && scatterlines!(ax1, r.X, r.out[m]; color=_o1color(m), markersize=4, linewidth=1, label=m)
        end
        axislegend(ax1; position=:rb, nbanks=2, labelsize=10)
        ax2 = Axis(fig[2, 1]; xlabel="x", ylabel="|error| (log)", yscale=log10)
        for m in methods
            lines!(ax2, r.X, abs.(r.out[m] .- r.q) .+ 1e-16; color=_o1color(m),
                   linestyle=occursin("lumped", m) ? :dash : :solid, label=@sprintf("%s (%.1e)", m, r.err[m]))
        end
        axislegend(ax2; position=:rb, nbanks=2, labelsize=9)
        appSave(fig, dir, @sprintf("%s_G%.1f", name, r.G))
        writedlm(joinpath(dir, @sprintf("%s_G%.1f.csv", name, r.G)),
                 vcat(permutedims(vcat(["x", "exact"], methods)), hcat(r.X, r.q, [r.out[m] for m in methods]...)), ',')
        global outputs1DLast = r
    elseif task == "convergence"
        c = outputs1DConvergence(name; L=Float64(A.L))
        dir = appOutputDir(@__DIR__, "convergence")
        fig = Figure(size=(950, 700))
        ax = Axis(fig[1, 1]; xlabel="points per wavelength G", ylabel="relative L2 error of q",
                  xscale=log10, yscale=log10, title="q = $(OUTPUT_REGISTRY[name].expr): output recipes")
        for m in c.methods
            scatterlines!(ax, c.Gs, c.err[m]; color=_o1color(m), label=m, linestyle=occursin("lumped", m) ? :dash : :solid)
        end
        for p in (2, 4, 6)
            lines!(ax, c.Gs, 0.3 .* (c.Gs[1] ./ c.Gs) .^ p; color=:gray85, linestyle=:dot)
        end
        Legend(fig[1, 2], ax; labelsize=10)
        appSave(fig, dir, "convergence_$(name)")
        writedlm(joinpath(dir, "convergence_$(name).csv"),
                 vcat(permutedims(vcat(["G"], c.methods)), hcat(c.Gs, [c.err[m] for m in c.methods]...)), ',')
        global outputs1DConvergenceLast = c
    else
        error("unknown task '$task' (compare, convergence, custom)")
    end
end

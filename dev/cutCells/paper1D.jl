# dev/cutCells/paper1D.jl -- 1-D free-end benchmark for the paper
#
# Elastic bar (famous equation "1DsismoTime"), free end at x_s = η Δx beyond the last
# grid node, point force (Gaussian weights) 2.5 wavelengths before the end, exact
# solution = direct wave + even image. Four free-surface treatments, all with the
# same weak-form rows, void extension and recipes:
#     OPT3 / OPT5  ×  zigzag (end snapped to the last node) / cut cell (exact end)
# Everything compared RAW (no normalisation, no amplitude fit).
#
#   (a) waveforms at one receiver for η = 0.5 (16 points per wavelength)
#   (b) reflection timing error vs η          (zigzag: exactly -2ηΔx/c)
#   (c) raw relative L2 error vs η
#   (d) raw relative L2 error vs points per wavelength (η = 0.37): convergence order
#
# Run:  julia --project=. dev/cutCells/paper1D.jl
# Output: dev/cutCells/outputs/paper1D/ (figure PNG + PDF, table CSV)

if !@isdefined(flexopt_root)
    import Pkg
    global flexopt_root = abspath(joinpath(@__DIR__, "..", ".."))
    Pkg.activate(flexopt_root)
end
include(joinpath(flexopt_root, "dev", "common", "devBootstrap.jl"))
include(joinpath(@__DIR__, "env.jl"))
using CairoMakie, Printf, LinearAlgebra, DelimitedFiles
CairoMakie.activate!()

base1D = Dict{String,Any}(
    "famousEquationType" => "1DsismoTime",
    "Δ" => (1.0, 1.0), "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3, "supplementaryOrder" => 2,
    "variationalForm" => DEFAULT_FORM[],
    "taylor_inverse_mode" => :weak_operator_optimized,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "recipe_backend" => backend,
)
schemes = Dict{String,Any}()
for s in ("OPT3", "OPT5")
    p = schemeParameters(base1D, s)
    schemes[s] = (params=p, recipe=cachedRecipe(p, "elastic1D_$s"), cut=flexOPT.CutCellScheme(p),
        reach=s == "OPT5" ? 2 : 1, extension=s == "OPT5" ? :cubic : :linear)
end

"One 1-D run; returns times, the receiver trace and the exact traces."
function freeEnd1DRun(s, geometry, η; ppw=16, receiverWavelengths=1.25, dx=1.0, c=1.0, rho=1.0, cfl=0.1)
    sc = schemes[s]
    λ = ppw * dx
    f0 = c / λ
    t0 = 1.5 / f0
    sourceDistance = 2.5λ
    receiverDistance = receiverWavelengths * λ
    xs = η * dx
    x = collect(-(sourceDistance + 4λ):dx:6dx)
    nx = length(x)
    duration = t0 + (sourceDistance + receiverDistance) / c + 1.2 / f0
    dt = cfl * dx / c
    ricker(t) = (a = π * f0 * (t - t0); (1 - 2a^2) * exp(-a^2))
    S(t) = (a = π * f0 * (t - t0); (t - t0) * exp(-a^2))
    ls = flexOPT.levelSetFromFunction(r -> r[1] - xs, (x,))
    geom = flexOPT.CutCellGeometry(ls; maxDegree=sc.cut.maxDegree, zigzag=(geometry == "zigzag"))
    solid = flexOPT.sideMask(geom.levelSet)
    cut = flexOPT.cutNodes(geom)
    aggregation = flexOPT.sliverAggregation(geom; extension=sc.extension, reach=sc.reach)
    extended = Set(a.node for a in aggregation)
    farVoid = [I for I in CartesianIndices(geom.levelSet.φ) if !solid[I] && !(I in cut) && !(I in extended)]
    models = [fill(1.0, nx), fill(c^2 * (dt / dx)^2, nx)]
    tm = sc.recipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching
    family = (models=models, modelPoints=getModelPoints(models[1], 3, tm), Δ=(1.0, 1.0),
        modelName="paper1D_$(s)_ppw$(ppw)")
    numerical = numericalOperatorConstruction(Dict{String,Any}("optRec" => sc.recipe, "modelFam" => family,
        "absorbingBoundaries" => nothing, "maskedRegionInSpace" => nothing, "boundaryConditions" => nothing,
        "representation" => "matrixfree"))["numOperators"]
    volume = prepareLinearSystem(numerical)
    rows = flexOPT.cutCellRows(sc.cut, geom, cut, models, volume)
    system = flexOPT.pinRows(flexOPT.applyAggregation(flexOPT.installRows(volume, rows), aggregation), farVoid)
    xsrc = xs - sourceDistance
    w = [exp(-(xi - xsrc)^2 / 2) for xi in x]
    w[abs.(x .- xsrc) .> 3] .= 0.0
    w ./= sum(w)
    steps = ceil(Int, duration / dt)
    nPast = system.timePointsUsedForOneStep - 1
    st = ((1-nPast):steps) .* dt
    source = zeros(system.NforcePoints, system.NForceField, length(st))
    for i in eachindex(w)
        w[i] == 0 && continue
        source[i, 1, :] .= w[i] * dt^2 / (rho * dx) .* ricker.(st)
    end
    prop = propagateLinearSystem(system, steps, dt; sourceFull=source, output_stride=max(1, round(Int, 0.02λ / dt)),
        blowup_limit=1e8, solver_name="paper1D", scheme=:direct)
    prop.stopped_early && error("blow-up: $s $geometry η=$η ppw=$ppw")
    ir = argmin(abs.(x .- (xs - receiverDistance)))
    xr = x[ir]
    times = Float64.(prop.times)
    numeric = Float64.(prop.history[ir, 1, :])
    direct(t) = sum(w[j] * S(t - abs(xr - x[j]) / c) for j in eachindex(x) if w[j] != 0) / (2rho * c)
    image(t) = sum(w[j] * S(t - abs(xr - (2xs - x[j])) / c) for j in eachindex(x) if w[j] != 0) / (2rho * c)
    exact = [direct(t) + image(t) for t in times]
    # reflection timing by cross-correlation of (numeric - exact direct) with the exact image
    tR = t0 + (2xs - xsrc - xr) / c
    tw = range(tR - 0.8 / f0, tR + 0.8 / f0; length=401)
    sampleN(t) = (k = clamp(searchsortedlast(times, t), 1, length(times) - 1);
        a = clamp((t - times[k]) / (times[k+1] - times[k]), 0.0, 1.0); (1 - a) * numeric[k] + a * numeric[k+1])
    ref = [image(t) for t in tw]
    lags = -2.0:0.005:2.0
    scores = [sum((sampleN(t + l) - direct(t + l)) * e for (t, e) in zip(tw, ref)) for l in lags]
    lag = lags[argmax(scores)]
    reflected = [sampleN(t) - direct(t) for t in tw]
    (times=times, numeric=numeric, exact=exact, xr=xr, lag=lag, λ=λ,
        peakRatio=maximum(abs, reflected) / maximum(abs, ref),
        relL2=norm(numeric .- exact) / norm(exact))
end

variants = [("OPT3", "zigzag"), ("OPT5", "zigzag"), ("OPT3", "cut"), ("OPT5", "cut")]
labelOf(s, g) = g == "cut" ? "$s cut cell" : "$s zigzag"
colorOf = Dict(("OPT3", "zigzag") => :darkorange, ("OPT5", "zigzag") => :purple,
    ("OPT3", "cut") => :dodgerblue, ("OPT5", "cut") => :crimson)
styleOf = Dict("zigzag" => :dash, "cut" => :solid)
markerOf = Dict("OPT3" => :circle, "OPT5" => :rect)

etas = collect(0.0:0.05:0.95)
ppws = [8, 10, 12, 16, 20, 24, 32]
etaConv = 0.37
table = Any[["variant" "eta" "ppw" "lag" "peakRatio" "relL2"]]
sweep = Dict{Any,Any}()
for (s, g) in variants
    for η in etas
        r = freeEnd1DRun(s, g, η)
        sweep[(s, g, :eta, η)] = r
        push!(table, [labelOf(s, g) η 16 r.lag r.peakRatio r.relL2])
    end
    for p in ppws
        r = freeEnd1DRun(s, g, etaConv; ppw=p)
        sweep[(s, g, :ppw, p)] = r
        push!(table, [labelOf(s, g) etaConv p r.lag r.peakRatio r.relL2])
    end
    @printf("%-16s  η sweep: max |lag| %.3f, max |peak-1| %.4f, max relL2 %.3e;  ppw 16→32 relL2 %.3e → %.3e\n",
        labelOf(s, g), maximum(abs(sweep[(s, g, :eta, η)].lag) for η in etas),
        maximum(abs(sweep[(s, g, :eta, η)].peakRatio - 1) for η in etas),
        maximum(sweep[(s, g, :eta, η)].relL2 for η in etas),
        sweep[(s, g, :ppw, 16)].relL2, sweep[(s, g, :ppw, 32)].relL2)
end
println("\nobserved convergence order of the raw rel. L2 error (η = $etaConv, ppw 16 → 32):")
for (s, g) in variants
    e1, e2 = sweep[(s, g, :ppw, 16)].relL2, sweep[(s, g, :ppw, 32)].relL2
    @printf("  %-16s %.2f\n", labelOf(s, g), log2(e1 / e2))
end

outdir = joinpath(@__DIR__, "outputs", "paper1D")
mkpath(outdir)
writedlm(joinpath(outdir, "paper1D_table.csv"), vcat(table...), ',')

fig = Figure(size=(1200, 900), fontsize=15)
# (a) waveforms
axA = Axis(fig[1, 1:2]; title="(a) receiver $(round(Int, 1.25 * 16)) Δx before the end, free end η = 0.5 Δx beyond the last node, 16 points/λ",
    xlabel="time (Δx / c)", ylabel="u (raw)")
ref = sweep[("OPT5", "cut", :eta, 0.5)]
lines!(axA, ref.times, ref.exact; color=:black, linewidth=4, label="exact")
for (s, g) in variants
    r = sweep[(s, g, :eta, 0.5)]
    lines!(axA, r.times, r.numeric; color=colorOf[(s, g)], linestyle=styleOf[g], linewidth=2, label=labelOf(s, g))
end
axislegend(axA; position=:lt, nbanks=2)
# (b) lag vs η
axB = Axis(fig[2, 1]; title="(b) reflection timing error", xlabel="η (fraction of Δx)", ylabel="lag (Δx / c)")
lines!(axB, [0, 1], [0, -2]; color=:gray, linestyle=:dot, label="-2η (snapped end)")
for (s, g) in variants
    scatterlines!(axB, etas, [sweep[(s, g, :eta, η)].lag for η in etas]; color=colorOf[(s, g)],
        linestyle=styleOf[g], marker=markerOf[s], label=labelOf(s, g))
end
axislegend(axB; position=:lb)
# (c) error vs η
axC = Axis(fig[2, 2]; title="(c) raw relative L2 error", xlabel="η (fraction of Δx)", ylabel="‖u - u_exact‖ / ‖u_exact‖", yscale=log10)
for (s, g) in variants
    scatterlines!(axC, etas, [max(sweep[(s, g, :eta, η)].relL2, 1e-8) for η in etas]; color=colorOf[(s, g)],
        linestyle=styleOf[g], marker=markerOf[s])
end
# (d) convergence
axD = Axis(fig[1:2, 3]; title="(d) convergence, η = $etaConv", xlabel="points per wavelength",
    ylabel="raw relative L2 error", xscale=log10, yscale=log10, xticks=(ppws, string.(ppws)))
for (s, g) in variants
    scatterlines!(axD, ppws, [sweep[(s, g, :ppw, p)].relL2 for p in ppws]; color=colorOf[(s, g)],
        linestyle=styleOf[g], marker=markerOf[s], label=labelOf(s, g))
end
for (order, y0) in ((1, 0.3), (2, 0.05), (4, 3e-3))
    lines!(axD, [8, 32], [y0, y0 * (8 / 32)^order]; color=:gray, linestyle=:dot)
    text!(axD, 8.4, y0 * 1.1; text="order $order", align=(:left, :bottom), fontsize=12, color=:gray30)
end
colsize!(fig.layout, 3, Relative(0.34))
save(joinpath(outdir, "paper1D_freeEnd.png"), fig; px_per_unit=2)
save(joinpath(outdir, "paper1D_freeEnd.pdf"), fig)
println("figure -> ", joinpath(outdir, "paper1D_freeEnd.png"), " (and .pdf); table -> paper1D_table.csv")

# dev/cutCells/tiltedPlane.jl -- step `tiltedPlane` of dev/cutCells/run.jl
#
# The validated flat point-force benchmark (dev/freeSurfacePointSource/run.jl:
# homogeneous half-space vp 6000, vs 3464, ρ 2700; Gaussian vertical force at
# 20 km depth, Ricker 0.65 Hz; dx 500 m; 9 s) rotated by θ on the SAME
# Cartesian grid: the free surface becomes the plane n̂·r = 0,
# n̂ = (-sin θ, cos θ), cut by the grid; source, force direction (along n̂)
# and receivers rotate with it. The exact answer is the flat solution
# rotated, so case tilt0 (same code, surface on grid nodes) is the
# reference; tilt10/20/30 measure what the cut-cell surface adds. Tilt0 is
# itself compared with the cached SPECFEM2D flat run when it is on disk.
#
# Receivers (flat frame): x' = -30, -20, -10, 0, 10 km at depth 10 km
# ("deep", the SPECFEM2D stations) and 1 km ("shallow", Rayleigh wave);
# components rotated back to (tangential u_t, normal u_n).
# No absorbing layer: every comparison stops before the first return from
# an artificial grid edge (computed for both grids, printed).
#
# Cached per case in data/dev/cutCells (prefix tiltedPlane); --redo=tilt
# recomputes. ENV FLEXOPT_TILT_DX changes Δx (500 m default).

include(joinpath(@__DIR__, "env.jl"))
using LinearAlgebra, StaticArrays, Statistics

const TP = (vp=6000.0, vs=3464.0, rho=2700.0, halfWidth=60e3, voidCapCells=3, duration=9.0,
    sourceFlat=(-10e3, -20e3), receiversX=(-30e3, -20e3, -10e3, 0.0, 10e3),
    lines=(deep=-10e3, shallow=-1e3), force=1.0e10, sigma=500.0, cutoff=2.0,
    outputSampling=0.10, referenceDx=500.0, pointsPerSWavelength=16, rickerSharpness=1.5,
    temporalRefinement=2)
TILT_CODE_VERSION = "v1-2026-11-24"
dxTilt = parse(Float64, get(ENV, "FLEXOPT_TILT_DX", "500.0"))
fTilt = TP.rickerSharpness * TP.vs / (TP.pointsPerSWavelength * TP.referenceDx)
delayTilt = 1.5 / fTilt
dtTilt = 0.20 * dxTilt / (sqrt(2) * TP.vp * TP.temporalRefinement)
rickerTilt(t) = (a = π * fTilt * (t - delayTilt); (1 - 2a^2) * exp(-a^2))
thetaOfCase(case) = parse(Float64, replace(case, "tilt" => ""))

optParameters2D = Dict{String,Any}(
    "famousEquationType" => "2DsismoTimeIsoHeteroSingleForce",
    "Δ" => (1.0, 1.0, 1.0), "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3, "supplementaryOrder" => 2,
    "variationalForm" => DEFAULT_FORM[],
    "taylor_inverse_mode" => :weak_operator_optimized,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "recipe_backend" => backend,
)
recipe2D = cachedRecipe(optParameters2D, "elastic2D_OPT3")
@time scheme2D = flexOPT.CutCellScheme(optParameters2D)
@show scheme2D.maxDegree scheme2D.pointsShape

function bilinearTrace(history, x, z, field, p)
    dx = x[2] - x[1]
    fi = (p[1] - x[1]) / dx
    fk = (p[2] - z[1]) / dx
    i, k = floor(Int, fi) + 1, floor(Int, fk) + 1
    a, b = fi - (i - 1), fk - (k - 1)
    values = (1 - a) * (1 - b) .* history[i, k, field, :] .+ a * (1 - b) .* history[i+1, k, field, :] .+
        (1 - a) * b .* history[i, k+1, field, :] .+ a * b .* history[i+1, k+1, field, :]
    nodes = [CartesianIndex(i, k), CartesianIndex(i + 1, k), CartesianIndex(i, k + 1), CartesianIndex(i + 1, k + 1)]
    return Float64.(values), nodes
end

function runTiltedPlane(θdeg)
    θ = deg2rad(θdeg)
    t̂ = SVector(cos(θ), sin(θ))
    n̂ = SVector(-sin(θ), cos(θ))
    toGrid(xp, zp) = xp * t̂ + zp * n̂
    dx, dt = dxTilt, dtTilt
    sx, sz = TP.sourceFlat
    W = TP.halfWidth
    corners = [toGrid(a, b) for a in (sx - W, sx + W), b in (sz - W, TP.voidCapCells * dx)]
    lo(d) = floor(minimum(c[d] for c in corners) / dx) * dx
    hi(d) = ceil(maximum(c[d] for c in corners) / dx) * dx
    x = collect(lo(1):dx:hi(1))
    z = collect(lo(2):dx:hi(2))
    nx, nz = length(x), length(z)
    println("  grid $(nx) × $(nz), Δx = $dx m, Δt = $(round(dt; sigdigits=4)) s")

    ls = flexOPT.levelSetFromFunction(r -> dot(n̂, r), (x, z))
    geom = flexOPT.CutCellGeometry(ls; maxDegree=scheme2D.maxDegree)
    solid = flexOPT.sideMask(ls)
    cut = flexOPT.cutNodes(geom)
    cutSet = Set(cut)
    aggregation = flexOPT.sliverAggregation(geom; extension=CUT_EXTENSION)
    extended = Set(a.node for a in aggregation)
    farVoid = [I for I in CartesianIndices(ls.φ) if !solid[I] && !(I in cutSet) && !(I in extended)]
    nMerged = count(a -> a.aggregate, aggregation)
    println("  $(length(cut)) cut nodes ($(count(I -> solid[I], cut)) solid), $(nMerged) void nodes merged, ",
        "$(length(aggregation) - nMerged) extended, $(length(farVoid)) pinned")

    # material: homogeneous, scaled as the volume benchmark, extended over the void
    rhoS = fill(1.0, nx, nz)
    lambdaS = fill((TP.vp^2 - 2TP.vs^2) * (dt / dx)^2, nx, nz)
    muS = fill(TP.vs^2 * (dt / dx)^2, nx, nz)
    models = [rhoS, lambdaS, muS]
    timeMarching = recipe2D["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching
    family = (models=models, modelPoints=getModelPoints(models[1], 3, timeMarching),
        Δ=(1.0, 1.0, 1.0), modelName="cutcells_tilt$(round(Int, θdeg))_dx$(round(Int, dx))")
    numerical = numericalOperatorConstruction(Dict{String,Any}(
        "optRec" => recipe2D, "modelFam" => family, "absorbingBoundaries" => nothing,
        "maskedRegionInSpace" => nothing, "boundaryConditions" => nothing,
        "representation" => "matrixfree"))["numOperators"]
    volume = prepareLinearSystem(numerical)
    volume.spaceShape == (nx, nz) || error("operator grid $(volume.spaceShape) differs from ($nx, $nz)")

    src = toGrid(sx, sz)
    probe = CartesianIndex(argmin(abs.(x .- src[1])), argmin(abs.(z .- src[2])))
    selfCheck = flexOPT.cutCellSelfCheck(scheme2D, volume, models, probe)
    @show selfCheck
    maximum(values(selfCheck)) < selfCheckTolerance(recipe2D) || error("cut-cell rows of an uncut node differ from the volume rows: $selfCheck")
    rowsTime = @elapsed rows = flexOPT.cutCellRows(scheme2D, geom, cut, models, volume)
    println("  cut-cell rows: $(round(rowsTime; digits=1)) s")
    system = flexOPT.pinRows(flexOPT.applyAggregation(flexOPT.installRows(volume, rows), aggregation), farVoid)

    # source: Gaussian weights around the rotated source point, force along n̂
    specs = [(I, exp(-sum(abs2, SVector(x[I[1]], z[I[2]]) - src) / (2TP.sigma^2)))
             for I in CartesianIndices((nx, nz)) if norm(SVector(x[I[1]], z[I[2]]) - src) <= TP.cutoff * TP.sigma]
    total = sum(last, specs)
    all(solid[first(s)] for s in specs) || error("source touches the void")
    steps = ceil(Int, TP.duration / dt)
    nPast = system.timePointsUsedForOneStep - 1
    sourceTimes = ((1-nPast):steps) .* dt
    wavelet = rickerTilt.(sourceTimes)
    scale = TP.force * dt^2 / (TP.rho * dx^2)
    source = zeros(system.NforcePoints, system.NForceField, length(sourceTimes))
    linear = LinearIndices((nx, nz))
    for (I, w) in specs
        source[linear[I], 1, :] .+= (w / total) * scale * n̂[1] .* wavelet
        source[linear[I], 2, :] .+= (w / total) * scale * n̂[2] .* wavelet
    end
    stride = max(1, round(Int, TP.outputSampling / dt))
    marchTime = @elapsed propagation = propagateLinearSystem(system, steps, dt; sourceFull=source,
        output_stride=stride, blowup_limit=1e8, solver_name="OPT3 cut-cell tilt $(θdeg)°", scheme=:direct)
    propagation.stopped_early && error("propagation stopped early (blow-up) at θ = $θdeg°")
    println("  propagation: $(round(marchTime; digits=1)) s")
    history = propagation.history                 # (nx, nz, 2, nt)

    # receivers, rotated components, first return from an artificial edge
    edges = [SVector(x[I[1]], z[I[2]]) for I in CartesianIndices((nx, nz))
             if solid[I] && (I[1] in (1, nx) || I[2] in (1, nz))]
    start = delayTilt - 1 / fTilt
    traces = Dict{String,Any}()
    for (line, zp) in pairs(TP.lines), (ir, xp) in enumerate(TP.receiversX)
        p = toGrid(xp, zp)
        ux, nodesUsed = bilinearTrace(history, x, z, 1, p)
        uz, _ = bilinearTrace(history, x, z, 2, p)
        all(solid[I] for I in nodesUsed) || error("receiver $(line) $(xp) interpolates from the void")
        tReturn = start + minimum(norm(src - e) + norm(e - p) for e in edges) / TP.vp
        traces["$(line)_$(ir)"] = (x_km=xp / 1e3, line=String(line),
            ut=ux .* cos(θ) .+ uz .* sin(θ), un=-ux .* sin(θ) .+ uz .* cos(θ), tReturn=tReturn)
    end
    Dict{String,Any}("theta" => θdeg, "times" => Float64.(propagation.times), "traces" => traces,
        "selfCheck" => selfCheck, "nCut" => length(cut), "nMerged" => nMerged,
        "nExtended" => length(aggregation) - nMerged, "rowsTime" => rowsTime, "marchTime" => marchTime)
end

tiltCases = filter(c -> startswith(c, "tilt"), DEV.cases)
"tilt0" in tiltCases || pushfirst!(tiltCases, "tilt0")
resultsTilt = Dict{String,Any}()
for case in tiltCases
    θdeg = thetaOfCase(case)
    println("\n== tilted plane, θ = $(θdeg)°")
    params = Dict{String,Any}("stage" => "tiltedPlane", "theta" => θdeg, "dx" => dxTilt,
        "halfWidth" => TP.halfWidth, "duration" => TP.duration, "extension" => String(CUT_EXTENSION),
        "codeVersion" => TILT_CODE_VERSION)
    devRedo!("tilt", params, "dev/cutCells", "tiltedPlane")
    resultsTilt[case] = myProduceOrLoad(cfg -> runTiltedPlane(θdeg), params, "dev/cutCells", "tiltedPlane")
end

# ------------------------------------------------------------- comparison
function windowMetrics(times, a, ref, tEnd)
    keep = times .<= tEnd
    ta, aa, rr = times[keep], a[keep], ref[keep]
    scale = max(maximum(abs, rr), eps())
    relL2 = norm(aa .- rr) / max(norm(rr), eps())
    sample(v, t) = (k = clamp(searchsortedlast(ta, t), 1, length(ta) - 1);
        w = clamp((t - ta[k]) / (ta[k+1] - ta[k]), 0.0, 1.0); (1 - w) * v[k] + w * v[k+1])
    lags = -0.3:0.005:0.3
    inner = filter(t -> t + 0.3 <= ta[end] && t - 0.3 >= ta[1], ta)
    scores = [sum(sample(aa, t + lag) * sample(rr, t) for t in inner) for lag in lags]
    (relL2=relL2, maxErr=maximum(abs, aa .- rr) / scale, lag=lags[argmax(scores)],
        corr=dot(aa .- mean(aa), rr .- mean(rr)) / max(norm(aa .- mean(aa)) * norm(rr .- mean(rr)), eps()))
end

reference = resultsTilt["tilt0"]
println("\n== cut-cell surface vs the flat reference (tilt0), window ends before any artificial-edge return")
@printf("  %-7s %-8s %-7s %-4s %-9s %-10s %-10s %-8s %-8s\n", "case", "line", "x' km", "comp", "t_end s",
    "rel. L2", "max err", "lag s", "corr")
summaryTilt = []
for case in tiltCases
    case == "tilt0" && continue
    r = resultsTilt[case]
    for key in sort(collect(keys(r["traces"])))
        tr, rf = r["traces"][key], reference["traces"][key]
        tEnd = min(TP.duration, tr.tReturn, rf.tReturn)
        for comp in (:un, :ut)
            m = windowMetrics(r["times"], getproperty(tr, comp), getproperty(rf, comp), tEnd)
            @printf("  %-7s %-8s %-7.0f %-4s %-9.2f %-10.3e %-10.3e %+-8.3f %-8.4f\n", case, tr.line, tr.x_km,
                String(comp), tEnd, m.relL2, m.maxErr, m.lag, m.corr)
            push!(summaryTilt, (case=case, theta=r["theta"], line=tr.line, x_km=tr.x_km, comp=comp, m...))
        end
    end
end

# tilt0 against the cached SPECFEM2D flat run (deep line, u_z = u_n)
specfemOutputTilt = joinpath(flexopt_root, "data", "specfem2d_benchmarks",
    "homogeneous_flat_free_surface_dx500m_nrec23", "OUTPUT_FILES")
if dxTilt == 500.0 && isdir(specfemOutputTilt)
    !@isdefined(specfemBenchmark) && include(joinpath(flexopt_root, "src", "specfemBenchmark.jl"))
    stations = ["S0016", "S0022", "S0017", "S0023", "S0018"]      # x = -30, -20, -10, 0, 10 km at z = -10 km
    files = specfemBenchmark.find_specfem2d_traces(specfemOutputTilt; component=:z, seismotype=2,
        network="FX", stations=stations)
    byStation = Dict(split(basename(f), '.')[2] => f for f in files)
    println("\n== tilt0 (cut-cell flat surface) vs SPECFEM2D, deep line, u_z")
    specfemTraces = Dict{Int,Any}()
    for (ir, name) in enumerate(stations)
        raw = specfemBenchmark.read_specfem2d_trace(byStation[name]; time_shift=1.2 / fTilt)
        disp = zeros(length(raw.values))
        disp[2:end] .= cumsum(((raw.values[1:end-1] .+ raw.values[2:end]) ./ 2) .* diff(raw.time))
        specfemTraces[ir] = (time=Float64.(raw.time), values=disp)
        mine = (time=reference["times"], values=reference["traces"]["deep_$(ir)"].un)
        m = specfemBenchmark.waveform_metrics(specfemTraces[ir], mine; samples=1001)
        @printf("  x' = %4.0f km   correlation %.4f   relative error (after amplitude fit) %.3e   amplitude %.4f\n",
            TP.receiversX[ir] / 1e3, m.correlation, m.relative_error, m.optimal_amplitude)
    end
end

# ---------------------------------------------------------------- figures
colorsTilt = Dict("tilt10" => :dodgerblue, "tilt20" => :darkorange, "tilt30" => :crimson, "tilt45" => :purple)
for line in ("deep", "shallow")
    fig = Figure(size=(1500, 230 * length(TP.receiversX) + 80))
    Label(fig[0, 1:2], "Tilted free surface (OPT3 cut-cell), $(line) receivers ($(round(Int, -TP.lines[Symbol(line)] / 1e3)) km below the surface): " *
        "rotated components vs θ = 0 (black)"; fontsize=14)
    for (ir, xp) in enumerate(TP.receiversX), (ic, comp) in enumerate((:un, :ut))
        key = "$(line)_$(ir)"
        ax = Axis(fig[ir, ic]; ylabel="x' = $(round(Int, xp / 1e3)) km", title=ir == 1 ? (comp == :un ? "u_n (normal)" : "u_t (tangential)") : "",
            xlabel=ir == length(TP.receiversX) ? "time (s)" : "")
        lines!(ax, reference["times"], getproperty(reference["traces"][key], comp); color=:black, linewidth=2.5, label="θ = 0°")
        for case in tiltCases
            case == "tilt0" && continue
            r = resultsTilt[case]
            lines!(ax, r["times"], getproperty(r["traces"][key], comp); color=get(colorsTilt, case, :gray),
                linestyle=:dash, label="θ = $(round(Int, r["theta"]))°")
            vlines!(ax, [min(TP.duration, r["traces"][key].tReturn, reference["traces"][key].tReturn)];
                color=get(colorsTilt, case, :gray), linestyle=:dot, linewidth=1)
        end
        ir == 1 && ic == 1 && axislegend(ax; position=:lt)
    end
    devSave(fig, "tilted_plane_$(line)")
end
println("\n(dotted lines: end of the comparison window, first return from an artificial grid edge)")

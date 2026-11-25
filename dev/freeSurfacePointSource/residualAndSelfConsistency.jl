# dev/freeSurfacePointSource/residualAndSelfConsistency.jl -- two reviewer-facing
# additions to the point-force/flat-free-surface, all-solvers comparison
# (dev/freeSurfacePointSource/runWithReferences.jl, job 045's output), without
# re-running any of the expensive FD3/OPT3/OPT5/SPECFEM2D propagation:
#
# 1. A residual panel under each existing waveform panel: for every solver,
#    reference - optimal_amplitude*candidate (the SAME amplitude-corrected
#    residual that src/specfemBenchmark.jl's waveform_metrics() scores),
#    plotted as a percentage of the reference trace's own RMS. This turns an
#    abstract "3% relative_error" into something visible on the same axes as
#    the waveform itself -- a reviewer can see directly that OPT5's 3% is a
#    near-flat, low-amplitude residual, not a hidden gross mismatch.
#
# 2. A SPECFEM2D self-consistency "noise floor" check: data/specfem2d_benchmarks
#    already has this exact physical case (homogeneous half-space, same
#    source/frequency) independently meshed at THREE grid spacings --
#    dx500m_nrec23 (this benchmark's reference), dx250m and dx375m (both
#    confirmed to already include the SAME 5 stations at z=-10km, just under
#    different station numbering -- see their DATA/STATIONS). Comparing
#    dx500m_nrec23 against dx250m/dx375m with the SAME waveform_metrics()
#    used for every solver gives a relative_error for "SPECFEM2D vs itself at
#    a different resolution" -- i.e. what a reviewer should expect even from
#    a textbook-perfect numerical match. No new simulation is needed; this is
#    pure post-processing of already-cached SPECFEM2D output.
#
# Run:
#   julia --project=. dev/freeSurfacePointSource/residualAndSelfConsistency.jl run
# REPL:
#   julia> include("dev/freeSurfacePointSource/residualAndSelfConsistency.jl")

if !@isdefined(flexopt_root)
    import Pkg
    function find_flexopt_root(start_dir=@__DIR__)
        candidates = String[]
        haskey(ENV, "FLEXOPT_ROOT") && push!(candidates, abspath(expanduser(ENV["FLEXOPT_ROOT"])))
        directory = abspath(start_dir)
        while true
            push!(candidates, directory)
            parent = dirname(directory)
            parent == directory && break
            directory = parent
        end
        for candidate in unique(candidates)
            isfile(joinpath(candidate, "src", "flexOPT.jl")) && return candidate
        end
        error("Cannot locate flexOPT; set ENV[\"FLEXOPT_ROOT\"]")
    end
    global flexopt_root = find_flexopt_root()
    Pkg.activate(flexopt_root)
end
include(joinpath(flexopt_root, "src", "specfemBenchmark.jl"))
using .specfemBenchmark
using CairoMakie, Statistics, LinearAlgebra, JLD2
CairoMakie.activate!(type="png")

# ----------------------------------------------------- load the cached run
tracesPath = joinpath(flexopt_root, "dev", "freeSurfacePointSource", "outputs",
    "run_with_references", "all_solvers_traces.jld2")
isfile(tracesPath) || error("Run dev/freeSurfacePointSource/runWithReferences.jl first -- $tracesPath not found")
cached = load(tracesPath)
receiverX = cached["receiverX"]
specfemWaveformsZ = cached["specfemWaveformsZ"]
sourceForce, sourceFrequency = cached["sourceForce"], cached["sourceFrequency"]

function solver_traces_for_station(station)
    base = (FD3=(time=cached["fdDisplacementZ"].time, values=cached["fdDisplacementZ"].values[:, station]),
        convFD3=(time=cached["convFDDisplacementZ"].time, values=cached["convFDDisplacementZ"].values[:, station]),
        OPT3=(time=cached["opt3DisplacementZ"].time, values=cached["opt3DisplacementZ"].values[:, station]),
        OPT5=(time=cached["opt5DisplacementZ"].time, values=cached["opt5DisplacementZ"].values[:, station]))
    get(cached, "haveSeismicWavePC", false) ?
        merge(base, (SeismicWavePC=cached["seismicWavePCDisplacementZ"][station],)) : base
end
haveSeismicWavePC = get(cached, "haveSeismicWavePC", false)

# ------------------------------------------- residual trace (same math as
# src/specfemBenchmark.jl's waveform_metrics, but returning the residual
# array itself rather than only scalar summaries -- duplicated locally
# rather than editing src/ without permission).
function _linear_sample_local(time, values, query)
    query <= first(time) && return first(values)
    query >= last(time) && return last(values)
    i = searchsortedlast(time, query)
    α = (query - time[i]) / (time[i + 1] - time[i])
    (1 - α) * values[i] + α * values[i + 1]
end
function residual_trace(reference, candidate; samples=2001)
    start_time = max(first(reference.time), first(candidate.time))
    end_time = min(last(reference.time), last(candidate.time))
    time = collect(range(start_time, end_time; length=samples))
    a = [_linear_sample_local(reference.time, reference.values, t) for t in time]
    b = [_linear_sample_local(candidate.time, candidate.values, t) for t in time]
    a .-= sum(a) / length(a)
    b .-= sum(b) / length(b)
    amplitude = sum(abs2, b) == 0 ? 0.0 : sum(a .* b) / sum(abs2, b)
    referenceRMS = sqrt(sum(abs2, a) / length(a))
    residualPercent = referenceRMS == 0 ? zero(a) : 100 .* (a .- amplitude .* b) ./ referenceRMS
    relativeError = referenceRMS == 0 ? NaN : sqrt(sum(abs2, a .- amplitude .* b) / length(a)) / referenceRMS
    (time=time, residualPercent=residualPercent, amplitude=amplitude, relativeError=relativeError)
end

# ------------------------------------------------- waveform + residual figure
solverColors = Dict(:FD3 => :steelblue, :convFD3 => :orange, :OPT3 => :seagreen,
    :OPT5 => :crimson, :SeismicWavePC => :purple)
combinedFigure = Figure(size=(1200, 330 * length(receiverX)))
for (row, station) in enumerate(eachindex(receiverX))
    waveAxis = Axis(combinedFigure[2row - 1, 1]; ylabel="x=$(receiverX[station]/1e3) km\nnormalized",
        xlabel="")
    residAxis = Axis(combinedFigure[2row, 1]; ylabel="residual (% of\nreference RMS)",
        xlabel=(station == lastindex(receiverX) ? "time (s)" : ""))
    reference = specfemWaveformsZ[station]
    traces = merge(solver_traces_for_station(station), (SPECFEM2D=reference,))
    for (name, trace) in pairs(traces)
        scale = max(maximum(abs, trace.values), eps(Float64))
        if name === :SPECFEM2D
            lines!(waveAxis, trace.time, trace.values ./ scale; label=string(name), color=:black, linewidth=2)
        else
            lines!(waveAxis, trace.time, trace.values ./ scale; label=string(name),
                color=get(solverColors, name, :gray),
                linestyle=(name === :SeismicWavePC ? :dot : :solid))
            r = residual_trace(reference, trace)
            lines!(residAxis, r.time, r.residualPercent; label=string(name),
                color=get(solverColors, name, :gray),
                linestyle=(name === :SeismicWavePC ? :dot : :solid))
        end
    end
    station == 1 && axislegend(waveAxis; position=:rb, nbanks=3)
end
Label(combinedFigure[0, 1],
    "Point force, flat free surface, u_z: waveform fit (top) and amplitude-corrected residual (bottom), " *
    "vs SPECFEM2D"; fontsize=14)
outdir = joinpath(flexopt_root, "dev", "freeSurfacePointSource", "outputs", "run_with_references")
mkpath(outdir)
save(joinpath(outdir, "waveform_fit_with_residuals.png"), combinedFigure)
println("saved figure -> ", joinpath(outdir, "waveform_fit_with_residuals.png"))

# --------------------------------------------- SPECFEM2D self-consistency
# Same 5 physical stations (x=-30,-20,-10,0,10 km @ z=-10km), read from two
# independently-meshed SPECFEM2D cases at different dx, compared against the
# SAME reference (dx500m_nrec23) this whole benchmark uses -- with the
# SAME waveform_metrics() scoring used for every solver above.
specfemTimeShift = 1.2 / sourceFrequency
function integrate_velocity_trace(trace)
    displacement = zeros(Float64, length(trace.values))
    displacement[2:end] .= cumsum(((trace.values[1:end-1] .+ trace.values[2:end]) ./ 2) .* diff(trace.time))
    (time=Float64.(trace.time), values=displacement)
end
function read_case_waveforms(caseDir, stationNames)
    output = joinpath(flexopt_root, "data", "specfem2d_benchmarks", caseDir, "OUTPUT_FILES")
    files = find_specfem2d_traces(output; component=:z, seismotype=2, network="FX", stations=stationNames)
    @assert length(files) == length(stationNames)
    byStation = Dict(split(basename(path), '.')[2] => path for path in files)
    ordered = [byStation[name] for name in stationNames]
    raw = [read_specfem2d_trace(file; time_shift=specfemTimeShift) for file in ordered]
    integrate_velocity_trace.(raw)
end
# x=-30,-20,-10,0,10 km @ z=-10km, in receiverX order, per each case's own
# DATA/STATIONS (confirmed by inspection before writing this script).
dx250mStations = ["S0007", "S0013", "S0008", "S0014", "S0009"]
dx375mStations = ["S0007", "S0013", "S0008", "S0014", "S0009"]
selfConsistencyResults = Dict{String,Any}()
for (label, caseDir, stationNames) in (
    ("dx250m", "homogeneous_flat_free_surface_dx250m", dx250mStations),
    ("dx375m", "homogeneous_flat_free_surface_dx375m", dx375mStations),
)
    candidateWaveforms = read_case_waveforms(caseDir, stationNames)
    perStation = [waveform_metrics(specfemWaveformsZ[s], candidateWaveforms[s]; samples=1001)
        for s in eachindex(receiverX)]
    summary = (mean_correlation=mean(m.correlation for m in perStation),
        rms_relative_error=sqrt(mean(abs2, (m.relative_error for m in perStation))),
        maximum_relative_error=maximum(m.relative_error for m in perStation))
    selfConsistencyResults[label] = summary
    println("SPECFEM2D self-consistency, dx500m_nrec23 (reference) vs ", label, ": ", summary)
end
display(selfConsistencyResults)
jldsave(joinpath(outdir, "specfem_self_consistency.jld2"); selfConsistencyResults)
println("saved -> ", joinpath(outdir, "specfem_self_consistency.jld2"))
println("DONE")

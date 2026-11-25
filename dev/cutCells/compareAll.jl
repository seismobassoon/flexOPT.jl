# dev/cutCells/compareAll.jl -- all six codes on the topography cases, RAW amplitudes
#
#   SPECFEM2D (reference, mesh follows the surface)      black
#   OPT3 cut-cell / OPT3 zigzag                          blue / orange dashed
#   OPT5 cut-cell / OPT5 zigzag                          red / purple dashed
#   SeismicWavePC (staircase, OpenSWPC-type FD)          green dotted
#
# Cases tilt0 (flat, all codes on the same problem), tilt10, tilt20, tilt30, hill.
# Every flexOPT/SPECFEM2D trace is raw (same force 1e10, same Gaussian source, same
# Ricker). SeismicWavePC has its own force normalisation and Ricker convention, so ONE
# global factor and ONE lag are fitted on tilt0 against SPECFEM2D (all traces at once) and
# applied unchanged to every case; the values are printed and written on the figures.
# Errors: rel. L2 against SPECFEM2D, both components together (|u| norm), mean over the
# five receivers of a line, window 0 ... min(9 s, first return from flexOPT's grid edge).
#
# Inputs (any missing one is skipped):
#   flexOPT caches  data/dev/cutCells/tiltedPlane_<hash>.jld2  (run.jl tiltedPlane, 4 variants)
#   SPECFEM2D       dev/cutCells/outputs/references/specfem_<case>.jld2        (specfemCases.jl)
#   SeismicWavePC   dev/cutCells/outputs/references/seismicwavepc_<case>.jld2  (seismicWavePCCases.jl)
# Outputs: dev/cutCells/outputs/compareAll/wiggles_<case>.png, errors.csv
#
# Run:  julia --project=. dev/cutCells/compareAll.jl [tilt10 tilt20 tilt30 hill]

if !@isdefined(flexopt_root)
    import Pkg
    global flexopt_root = abspath(joinpath(@__DIR__, "..", ".."))
    Pkg.activate(flexopt_root)
end
!@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
using .commonBatchs
include(joinpath(@__DIR__, "geometryCases.jl"))
using DrWatson, JLD2, CairoMakie, Printf, LinearAlgebra, Statistics
CairoMakie.activate!(type="png")

const CA = (dx=500.0, codeVersion="v1-2026-11-24")
flexVariants = [("OPT3", "cut"), ("OPT3", "zigzag"), ("OPT5", "cut"), ("OPT5", "zigzag")]
codeLabel = Dict(("OPT3", "cut") => "OPT3 cut-cell", ("OPT3", "zigzag") => "OPT3 zigzag",
    ("OPT5", "cut") => "OPT5 cut-cell", ("OPT5", "zigzag") => "OPT5 zigzag",
    "specfem" => "SPECFEM2D", "swpc" => "SeismicWavePC (staircase)")
codeStyle = Dict(("OPT3", "cut") => (color=:dodgerblue, linestyle=:solid, linewidth=1.5),
    ("OPT3", "zigzag") => (color=:darkorange, linestyle=:dash, linewidth=1.5),
    ("OPT5", "cut") => (color=:crimson, linestyle=:solid, linewidth=1.5),
    ("OPT5", "zigzag") => (color=:purple, linestyle=:dash, linewidth=1.5),
    "specfem" => (color=:black, linestyle=:solid, linewidth=3.0),
    "swpc" => (color=:forestgreen, linestyle=:dot, linewidth=1.8))

# ---------------------------------------------------------------- loading
function flexRun(case, scheme, geometry)
    params = Dict{String,Any}("stage" => "tiltedPlane", "theta" => caseTheta(case), "dx" => CA.dx,
        "halfWidth" => GC.halfWidth, "duration" => GC.duration,
        "extension" => scheme == "OPT5" ? "cubic" : "linear", "codeVersion" => CA.codeVersion)
    caseKind(case) == :hill && (params["case"] = "hill")
    scheme == "OPT3" || (params["scheme"] = scheme)
    geometry == "cut" || (params["geometry"] = geometry)
    file = datadir("dev/cutCells", "tiltedPlane_$(commonBatchs.hash_parameters(params)).jld2")
    isfile(file) ? load(file) : nothing
end
function referenceRun(code, case)
    file = joinpath(@__DIR__, "outputs", "references", "$(code)_$(case).jld2")
    isfile(file) ? load(file) : nothing
end

function sampleAt(times, v, t)
    (t < times[1] || t > times[end]) && return 0.0
    k = clamp(searchsortedlast(times, t), 1, length(times) - 1)
    w = clamp((t - times[k]) / (times[k+1] - times[k]), 0.0, 1.0)
    (1 - w) * v[k] + w * v[k+1]
end
resample(times, v, tt; factor=1.0, lag=0.0) = [factor * sampleAt(times, v, t - lag) for t in tt]

caseList = isempty(ARGS) ? ["tilt10", "tilt20", "tilt30", "hill"] : ARGS
allCases = unique(vcat(["tilt0"], caseList))
data = Dict{String,Any}()
for case in allCases
    d = Dict{Any,Any}()
    for v in flexVariants
        r = flexRun(case, v...)
        isnothing(r) || (d[v] = r)
    end
    for code in ("specfem", "seismicwavepc")
        r = referenceRun(code, case)
        isnothing(r) || (d[code == "specfem" ? "specfem" : "swpc"] = r)
    end
    data[case] = d
    println("  $case: ", join([string(codeLabel[k]) for k in keys(d)], ", "))
end
lineNames = ("deep", "shallow")
nRec = length(GC.receiversX)

# window per case: before the first artificial-edge return of any flexOPT run (absorbing in the others)
function windowEnd(case, key)
    d = data[case]
    t = GC.duration
    for v in flexVariants
        haskey(d, v) && (t = min(t, d[v]["traces"][key].tReturn))
    end
    t
end

# ------------------------------------------- SeismicWavePC global calibration on tilt0
swpcFactor, swpcLag = NaN, NaN
if haskey(data["tilt0"], "specfem") && haskey(data["tilt0"], "swpc")
    ref, sw = data["tilt0"]["specfem"], data["tilt0"]["swpc"]
    tt = collect(0.0:0.01:GC.duration)
    best = Ref((Inf, 0.0, 0.0))
    for lag in -0.6:0.005:0.6
        num, den, rr = 0.0, 0.0, 0.0
        for line in lineNames, ir in 1:nRec, comp in (:un, :ut)
            key = "$(line)_$ir"
            keep = tt .<= windowEnd("tilt0", key)
            a = resample(ref["times"], getproperty(ref["traces"][key], comp), tt[keep])
            b = resample(sw["times"], getproperty(sw["traces"][key], comp), tt[keep]; lag=lag)
            num += dot(a, b); den += dot(b, b); rr += dot(a, a)
        end
        f = num / den
        misfit = rr - 2f * num + f^2 * den
        misfit < best[][1] && (best[] = (misfit, f, lag))
    end
    global swpcFactor, swpcLag = best[][2], best[][3]
    @printf("\nSeismicWavePC global calibration on tilt0 vs SPECFEM2D: factor %.4e, lag %+.3f s (applied to every case)\n",
        swpcFactor, swpcLag)
else
    println("\nSeismicWavePC calibration skipped (tilt0 SPECFEM2D or SeismicWavePC missing): SeismicWavePC not plotted")
end

"Trace of `code` at `key`, component `comp`, on times tt (raw; SeismicWavePC calibrated)."
function codeTrace(case, code, key, comp, tt)
    r = data[case][code]
    v = getproperty(r["traces"][key], comp)
    code == "swpc" ? resample(r["times"], v, tt; factor=swpcFactor, lag=swpcLag) : resample(r["times"], v, tt)
end
codesOf(case) = [c for c in vcat(["specfem"], flexVariants, ["swpc"]) if haskey(data[case], c) &&
    !(c == "swpc" && isnan(swpcFactor))]

# ---------------------------------------------------------------- error table
outDir = joinpath(@__DIR__, "outputs", "compareAll")
mkpath(outDir)
rows = String["case,line,code,relL2_vs_specfem_mean,relL2_vs_specfem_max"]
println("\n== raw rel. L2 vs SPECFEM2D (both components, mean / max over the 5 receivers of a line)")
@printf("  %-7s %-8s %-28s %8s %8s\n", "case", "line", "code", "mean", "max")
for case in allCases
    haskey(data[case], "specfem") || continue
    for line in lineNames, code in codesOf(case)
        code == "specfem" && continue
        errs = map(1:nRec) do ir
            key = "$(line)_$ir"
            tt = collect(0.0:0.01:windowEnd(case, key))
            dn = codeTrace(case, code, key, :un, tt) .- codeTrace(case, "specfem", key, :un, tt)
            dt = codeTrace(case, code, key, :ut, tt) .- codeTrace(case, "specfem", key, :ut, tt)
            rn, rt = codeTrace(case, "specfem", key, :un, tt), codeTrace(case, "specfem", key, :ut, tt)
            sqrt(sum(abs2, dn) + sum(abs2, dt)) / sqrt(sum(abs2, rn) + sum(abs2, rt))
        end
        @printf("  %-7s %-8s %-28s %8.3f %8.3f\n", case, line, codeLabel[code], mean(errs), maximum(errs))
        push!(rows, @sprintf("%s,%s,%s,%.5f,%.5f", case, line, codeLabel[code], mean(errs), maximum(errs)))
    end
end
write(joinpath(outDir, "errors.csv"), join(rows, "\n") * "\n")
println("table -> ", joinpath(outDir, "errors.csv"))

# ---------------------------------------------------------------- wiggle plots
# one figure per case: rows = receiver lines, columns = (u_n, u_t); traces offset by the receiver
# position; ONE scale per panel row (from SPECFEM2D's max over both components), identical for all
# codes, so amplitudes stay raw and directly comparable
for case in allCases
    codes = codesOf(case)
    isempty(codes) && continue
    tt = collect(0.0:0.01:GC.duration)
    spacing = 10.0                                  # km between receivers
    title = caseKind(case) == :hill ? "Gaussian hill (H = 3 km, w = 5 km)" :
        "Tilted free surface θ = $(Int(caseTheta(case)))°"
    fig = Figure(size=(1500, 1250))
    Label(fig[0, 1:2], title * " -- raw displacement, same scale for every code" *
        (isnan(swpcFactor) ? "" : @sprintf("; SeismicWavePC × %.3e, lag %+.3f s (global, from tilt0)", swpcFactor, swpcLag));
        fontsize=15)
    for (il, line) in enumerate(lineNames)
        scaleCode = haskey(data[case], "specfem") ? "specfem" : first(codes)
        peak = maximum(maximum(abs, codeTrace(case, scaleCode, "$(line)_$ir", comp, tt))
                       for ir in 1:nRec, comp in (:un, :ut))
        gain = 0.45 * spacing / max(peak, eps())
        for (ic, comp) in enumerate((:un, :ut))
            ax = Axis(fig[il, ic]; xlabel="time (s)",
                ylabel=caseKind(case) == :hill ? "receiver x (km)" : "receiver x' (km, along the surface)",
                title="$(line) line, " * (comp == :un ? "u_n (normal)" : "u_t (tangential)") *
                      (caseKind(case) == :hill ? " = " * (comp == :un ? "u_z" : "u_x") : ""),
                yticks=collect(GC.receiversX ./ 1e3))
            for ir in 1:nRec
                key = "$(line)_$ir"
                y0 = GC.receiversX[ir] / 1e3
                for code in codes
                    st = codeStyle[code]
                    lines!(ax, tt, y0 .+ gain .* codeTrace(case, code, key, comp, tt);
                        color=st.color, linestyle=st.linestyle, linewidth=st.linewidth,
                        label=ir == 1 ? codeLabel[code] : nothing)
                end
                tEnd = windowEnd(case, key)
                tEnd < GC.duration && scatter!(ax, [tEnd], [y0]; marker=:vline, markersize=14, color=:gray40)
            end
            il == 1 && ic == 1 && axislegend(ax; position=:lb, framevisible=true, labelsize=11)
        end
        Label(fig[il, 3], @sprintf("scale:\n%.2e m\nper 4.5 km", peak); rotation=0, fontsize=11, tellheight=false)
    end
    Label(fig[3, 1:2], "grey tick: first return from flexOPT's (non-absorbing) grid edge -- errors are computed before it";
        fontsize=11)
    path = joinpath(outDir, "wiggles_$(case).png")
    save(path, fig)
    println("figure -> ", path)
end

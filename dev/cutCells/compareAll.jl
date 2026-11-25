# dev/cutCells/compareAll.jl -- all codes on the topography cases, RAW amplitudes
#
# REFERENCE ("exact", black)
#   tilt cases: the FLAT SPECFEM2D solution (case tilt0, surface/mesh aligned with the grid),
#     taken at the same receiver: the tilted problem is the flat one rotated (source, force,
#     receivers and surface rotate together), so its (u_t, u_n) are exactly the flat (u_x, u_z).
#     No tilted mesh is involved in the reference.
#   hill: SPECFEM2D on the hill mesh (no exact solution); the refined run (1 km elements,
#     FLEXOPT_SPECFEM_ELEMENT=1000) is used when present.
#   In both, the refined SPECFEM2D run replaces the 2 km one when it exists.
# OTHER CODES
#   SPECFEM2D on the tilted / hill mesh (2 km elements)  grey   -- tests SPECFEM's own meshing
#   OPT3 cut-cell / OPT3 zigzag                          blue / orange dashed
#   OPT5 cut-cell / OPT5 zigzag                          red / purple dashed
#   SeismicWavePC (staircase)                            green dotted
# At tilt0 cut-cell and zigzag are the same run (surface on grid nodes): the curves coincide.
#
# SeismicWavePC: ONE factor and ONE lag fitted on tilt0 against the reference, applied to every
# case (printed and written on the figures). Nothing else is rescaled.
# Errors: rel. L2 against the reference, both components (|u| norm), mean/max over the five
# receivers of a line, window 0 ... min(9 s, first return from flexOPT's grid edge).
# Figures: wiggles (u_n, u_t) + residuals (code - reference, x5) per case; the title lists every
# code that is MISSING (no cache found), and the run start prints the expected cache files.
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

const CA = (dx=500.0, codeVersion="v1-2026-11-24", residualGain=5.0)
flexVariants = [("OPT3", "cut"), ("OPT3", "zigzag"), ("OPT5", "cut"), ("OPT5", "zigzag")]
codeLabel = Dict(("OPT3", "cut") => "OPT3 cut-cell", ("OPT3", "zigzag") => "OPT3 zigzag",
    ("OPT5", "cut") => "OPT5 cut-cell", ("OPT5", "zigzag") => "OPT5 zigzag",
    "specfem" => "SPECFEM2D (case mesh, 2 km)", "swpc" => "SeismicWavePC (staircase)")
codeStyle = Dict(("OPT3", "cut") => (color=:dodgerblue, linestyle=:solid, linewidth=1.5),
    ("OPT3", "zigzag") => (color=:darkorange, linestyle=:dash, linewidth=1.5),
    ("OPT5", "cut") => (color=:crimson, linestyle=:solid, linewidth=1.5),
    ("OPT5", "zigzag") => (color=:purple, linestyle=:dash, linewidth=1.8),
    "reference" => (color=:black, linestyle=:solid, linewidth=3.0),
    "specfem" => (color=:gray55, linestyle=:solid, linewidth=2.0),
    "swpc" => (color=:forestgreen, linestyle=:dot, linewidth=1.8))

# ---------------------------------------------------------------- loading
function flexCacheFile(case, scheme, geometry)
    params = Dict{String,Any}("stage" => "tiltedPlane", "theta" => caseTheta(case), "dx" => CA.dx,
        "halfWidth" => GC.halfWidth, "duration" => GC.duration,
        "extension" => scheme == "OPT5" ? "cubic" : "linear", "codeVersion" => CA.codeVersion)
    caseKind(case) == :hill && (params["case"] = "hill")
    scheme == "OPT3" || (params["scheme"] = scheme)
    geometry == "cut" || (params["geometry"] = geometry)
    datadir("dev/cutCells", "tiltedPlane_$(commonBatchs.hash_parameters(params)).jld2")
end
referenceFile(code, case) = joinpath(@__DIR__, "outputs", "references", "$(code)_$(case).jld2")

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
missingCodes = Dict{String,Vector{String}}()
println("== inventory (found / MISSING, with the file looked for)")
for case in allCases
    d = Dict{Any,Any}()
    miss = String[]
    entries = vcat([(v, flexCacheFile(case, v...)) for v in flexVariants],
        [("specfem", referenceFile("specfem", case)), ("specfemFine", referenceFile("specfem", case * "_el1000m")),
         ("swpc", referenceFile("seismicwavepc", case))])
    for (key, file) in entries
        if isfile(file)
            d[key] = load(file)
        else
            key == "specfemFine" || push!(miss, get(codeLabel, key, string(key)))
        end
        @printf("  %-7s %-30s %-8s %s\n", case, get(codeLabel, key, key == "specfemFine" ? "SPECFEM2D refined (1 km)" : string(key)),
            isfile(file) ? "found" : "MISSING", relpath(file, flexopt_root))
    end
    data[case] = d
    missingCodes[case] = miss
end
lineNames = ("deep", "shallow")
nRec = length(GC.receiversX)

"Reference run and its trace key for (case, key)."
function referenceOf(case)
    if caseKind(case) == :tilt
        d = data["tilt0"]
        haskey(d, "specfemFine") && return (d["specfemFine"], "flat SPECFEM2D (1 km elements) rotated = exact")
        haskey(d, "specfem") && return (d["specfem"], "flat SPECFEM2D rotated = exact")
    else
        d = data[case]
        haskey(d, "specfemFine") && return (d["specfemFine"], "SPECFEM2D hill mesh, 1 km elements")
        haskey(d, "specfem") && return (d["specfem"], "SPECFEM2D hill mesh, 2 km elements")
    end
    (nothing, "none")
end

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
let (ref, _) = referenceOf("tilt0")
    if !isnothing(ref) && haskey(data["tilt0"], "swpc")
        sw = data["tilt0"]["swpc"]
        tt = collect(0.0:0.01:GC.duration)
        best = Ref((Inf, 0.0, 0.0))
        for lag in -1.5:0.005:1.5     # SeismicWavePC ricker peaks ~0.8 s earlier than gcRicker
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
        @printf("\nSeismicWavePC global calibration on tilt0: factor %.4e, lag %+.3f s (applied to every case)\n",
            swpcFactor, swpcLag)
    else
        println("\nSeismicWavePC calibration skipped (flat SPECFEM2D or SeismicWavePC tilt0 missing): not plotted")
    end
end

function codeTrace(case, code, key, comp, tt)
    if code == "reference"
        r = first(referenceOf(case))
        return resample(r["times"], getproperty(r["traces"][key], comp), tt)
    end
    r = data[case][code]
    v = getproperty(r["traces"][key], comp)
    code == "swpc" ? resample(r["times"], v, tt; factor=swpcFactor, lag=swpcLag) : resample(r["times"], v, tt)
end
"Codes compared with the reference (the case-mesh SPECFEM2D only when it is not the reference itself)."
function codesOf(case)
    ref, _ = referenceOf(case)
    out = Any[]
    haskey(data[case], "specfem") && data[case]["specfem"] !== ref && push!(out, "specfem")
    append!(out, [v for v in flexVariants if haskey(data[case], v)])
    haskey(data[case], "swpc") && !isnan(swpcFactor) && push!(out, "swpc")
    out
end

# ---------------------------------------------------------------- error table
outDir = joinpath(@__DIR__, "outputs", "compareAll")
mkpath(outDir)
rows = String["case,reference,line,code,relL2_mean,relL2_max"]
println("\n== raw rel. L2 against the reference (both components; mean / max over the 5 receivers)")
for case in allCases
    ref, refName = referenceOf(case)
    isnothing(ref) && (println("  $case: no reference"); continue)
    println("  $case -- reference: $refName")
    for line in lineNames, code in codesOf(case)
        errs = map(1:nRec) do ir
            key = "$(line)_$ir"
            tt = collect(0.0:0.01:windowEnd(case, key))
            rn, rt = codeTrace(case, "reference", key, :un, tt), codeTrace(case, "reference", key, :ut, tt)
            dn = codeTrace(case, code, key, :un, tt) .- rn
            dt = codeTrace(case, code, key, :ut, tt) .- rt
            sqrt(sum(abs2, dn) + sum(abs2, dt)) / sqrt(sum(abs2, rn) + sum(abs2, rt))
        end
        @printf("    %-8s %-30s %8.3f %8.3f\n", line, codeLabel[code], mean(errs), maximum(errs))
        push!(rows, @sprintf("%s,%s,%s,%s,%.5f,%.5f", case, refName, line, codeLabel[code], mean(errs), maximum(errs)))
    end
end
write(joinpath(outDir, "errors.csv"), join(rows, "\n") * "\n")
println("table -> ", joinpath(outDir, "errors.csv"))

# ---------------------------------------------------------------- wiggle plots
for case in allCases
    ref, refName = referenceOf(case)
    codes = codesOf(case)
    (isnothing(ref) && isempty(codes)) && continue
    tt = collect(0.0:0.01:GC.duration)
    spacing = 10.0
    title = caseKind(case) == :hill ? "Gaussian hill (H = 3 km, w = 5 km)" :
        "Tilted free surface θ = $(Int(caseTheta(case)))°"
    sub = "black: $refName.  Raw displacement, one scale per line for every code" *
        (isnan(swpcFactor) ? "" : @sprintf("; SeismicWavePC × %.3e, lag %+.3f s (global, tilt0)", swpcFactor, swpcLag)) *
        (case == "tilt0" ? ".  At θ = 0 cut-cell = zigzag (same run)." : "") *
        (isempty(missingCodes[case]) ? "" : "\nMISSING (no cache): " * join(missingCodes[case], ", "))
    fig = Figure(size=(2300, 1300))
    Label(fig[0, 1:4], title; fontsize=17, font=:bold)
    Label(fig[1, 1:4], sub; fontsize=12)
    plotCodes = isnothing(ref) ? codes : vcat(["reference"], codes)
    for (il, line) in enumerate(lineNames)
        scaleCode = first(plotCodes)
        peak = maximum(maximum(abs, codeTrace(case, scaleCode, "$(line)_$ir", comp, tt)) for ir in 1:nRec, comp in (:un, :ut))
        gain = 0.45 * spacing / max(peak, eps())
        for (ic, (comp, residual)) in enumerate(((:un, false), (:ut, false), (:un, true), (:ut, true)))
            (residual && isnothing(ref)) && continue
            cname = comp == :un ? "u_n" : "u_t"
            ax = Axis(fig[il+1, ic]; xlabel="time (s)", yticks=collect(GC.receiversX ./ 1e3),
                ylabel=ic == 1 ? (caseKind(case) == :hill ? "$line line, receiver x (km)" : "$line line, x' (km, along surface)") : "",
                title=residual ? "$line: $cname − reference (× $(Int(CA.residualGain)))" : "$line: $cname")
            for ir in 1:nRec
                key = "$(line)_$ir"
                y0 = GC.receiversX[ir] / 1e3
                refTrace = isnothing(ref) ? nothing : codeTrace(case, "reference", key, comp, tt)
                for code in plotCodes
                    residual && code == "reference" && continue
                    st = codeStyle[code]
                    v = codeTrace(case, code, key, comp, tt)
                    residual && (v = CA.residualGain .* (v .- refTrace))
                    lines!(ax, tt, y0 .+ gain .* v; color=st.color, linestyle=st.linestyle, linewidth=st.linewidth,
                        label=ir == 1 ? (code == "reference" ? "reference" : codeLabel[code]) : nothing)
                end
                tEnd = windowEnd(case, key)
                tEnd < GC.duration && scatter!(ax, [tEnd], [y0]; marker=:vline, markersize=16, color=:gray30)
            end
            il == 1 && ic == 1 && axislegend(ax; position=:lb, framevisible=true, labelsize=11)
        end
        Label(fig[il+1, 5], @sprintf("scale:\n%.2e m\n= 4.5 km", peak); fontsize=11, tellheight=false)
    end
    Label(fig[4, 1:4], "grey tick: first return from flexOPT's (non-absorbing) grid edge -- errors are computed before it"; fontsize=11)
    path = joinpath(outDir, "wiggles_$(case).png")
    save(path, fig)
    println("figure -> ", path)
end

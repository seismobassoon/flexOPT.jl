# dev/cutCells/compare3D.jl -- the 3-D cases, all codes, RAW amplitudes
#
# Reference (black): for every plane case the FLAT SPECFEM3D run (flat3D), taken at the same
# receiver -- the oblique problem is the flat one rotated, so its (u_1, u_2, u_n) are exactly
# the flat (u_x, u_y, u_z). SPECFEM3D on the oblique mesh is drawn in grey (tests its meshing).
# flexOPT: OPT3/OPT5 cut-cell (blue/red), zigzag (orange/purple dashed).
# Error: rel. L2 against the reference, the three components together, mean/max over the
# 9 receivers of a line, window 0-7 s (every code absorbs at its artificial boundaries).
# Figures: wiggles of u_n, u_1, u_2 and the residual (code − reference) ×5 of u_n and u_1;
# receivers 1-5 along x' (y' = 0), 6-9 along y' (x' = 0).
#
# Run:  julia --project=. dev/cutCells/compare3D.jl [oblique10 oblique20]   (add --small for the test size)

if !@isdefined(flexopt_root)
    import Pkg
    global flexopt_root = abspath(joinpath(@__DIR__, "..", ".."))
    Pkg.activate(flexopt_root)
end
!@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
using .commonBatchs
include(joinpath(@__DIR__, "geometryCases3D.jl"))
using DrWatson, JLD2, CairoMakie, Printf, Statistics
CairoMakie.activate!(type="png")

SMALL = "--small" in ARGS
cases = filter(a -> !startswith(a, "--"), ARGS)
isempty(cases) && (cases = ["flat3D", "oblique10", "oblique20"])
size3D = SMALL ? (halfWidth=5e3, below=5e3, padding=6, duration=2.0) :
    (halfWidth=GC3.halfWidth, below=GC3.below, padding=GC3.cerjanCells, duration=GC3.duration)
variants = [("OPT3", "cut"), ("OPT3", "zigzag"), ("OPT5", "cut"), ("OPT5", "zigzag")]
label = Dict(("OPT3", "cut") => "OPT3 cut-cell", ("OPT3", "zigzag") => "OPT3 zigzag",
    ("OPT5", "cut") => "OPT5 cut-cell", ("OPT5", "zigzag") => "OPT5 zigzag",
    "reference" => "reference: flat SPECFEM3D rotated (exact)", "specfem" => "SPECFEM3D, oblique mesh")
style = Dict(("OPT3", "cut") => (color=:dodgerblue, linestyle=:solid, linewidth=1.5),
    ("OPT3", "zigzag") => (color=:darkorange, linestyle=:dash, linewidth=1.5),
    ("OPT5", "cut") => (color=:crimson, linestyle=:solid, linewidth=1.5),
    ("OPT5", "zigzag") => (color=:purple, linestyle=:dash, linewidth=1.8),
    "reference" => (color=:black, linestyle=:solid, linewidth=3.0),
    "specfem" => (color=:gray55, linestyle=:solid, linewidth=2.0))

function flexFile(case, scheme, geometry)
    p = Dict{String,Any}("stage" => "run3D", "case" => case, "scheme" => scheme, "geometry" => geometry,
        "extension" => scheme == "OPT5" ? "cubic" : "linear", "dx" => GC3.dx,
        "halfWidth" => size3D.halfWidth, "below" => size3D.below, "padding" => size3D.padding,
        "duration" => size3D.duration, "codeVersion" => "v1-2026-11-26")
    datadir("dev/cutCells", "run3D_$(commonBatchs.hash_parameters(p)).jld2")
end
refFile(case) = joinpath(@__DIR__, "outputs", "references", "specfem3d_$(case).jld2")

function sampleAt(ts, v, t)
    (t < ts[1] || t > ts[end]) && return 0.0
    k = clamp(searchsortedlast(ts, t), 1, length(ts) - 1)
    w = clamp((t - ts[k]) / (ts[k+1] - ts[k]), 0.0, 1.0)
    (1 - w) * v[k] + w * v[k+1]
end

println("== inventory")
reference = isfile(refFile("flat3D")) ? load(refFile("flat3D")) : nothing
@printf("  %-10s %-38s %s\n", "flat3D", label["reference"], isnothing(reference) ? "MISSING" : "found")
data = Dict{String,Any}()
for case in cases
    d = Dict{Any,Any}()
    case != "flat3D" && isfile(refFile(case)) && (d["specfem"] = load(refFile(case)))
    case != "flat3D" && @printf("  %-10s %-38s %s\n", case, label["specfem"], haskey(d, "specfem") ? "found" : "MISSING")
    for v in variants
        f = flexFile(case, v...)
        isfile(f) && (d[v] = load(f))
        @printf("  %-10s %-38s %s  %s\n", case, label[v], isfile(f) ? "found  " : "MISSING", relpath(f, flexopt_root))
    end
    data[case] = d
end
isnothing(reference) && error("flat SPECFEM3D reference missing: run specfem3DCases.jl flat3D")

keysOf(line) = ["$(line)_$i" for i in 1:9]
tt = collect(0.0:0.01:size3D.duration)
trace(r, key, comp) = [sampleAt(r["times"], getproperty(r["traces"][key], comp), t) for t in tt]
codes(case) = vcat(haskey(data[case], "specfem") ? ["specfem"] : [], [v for v in variants if haskey(data[case], v)])

outDir = joinpath(@__DIR__, "outputs", "compare3D" * (SMALL ? "_small" : ""))
mkpath(outDir)
rows = ["case,line,code,relL2_mean,relL2_max"]
println("\n== raw rel. L2 against the flat SPECFEM3D reference (u_1, u_2, u_n together)")
for case in cases, line in ("deep", "shallow"), code in codes(case)
    errs = map(keysOf(line)) do key
        num = sum(sum(abs2, trace(data[case][code], key, c) .- trace(reference, key, c)) for c in (:u1, :u2, :un))
        den = sum(sum(abs2, trace(reference, key, c)) for c in (:u1, :u2, :un))
        sqrt(num / den)
    end
    @printf("  %-10s %-8s %-28s %7.3f %7.3f\n", case, line, label[code], mean(errs), maximum(errs))
    push!(rows, @sprintf("%s,%s,%s,%.5f,%.5f", case, line, label[code], mean(errs), maximum(errs)))
end
write(joinpath(outDir, "errors.csv"), join(rows, "\n") * "\n")

for case in cases
    cs = codes(case)
    fig = Figure(size=(2600, 1400))
    Label(fig[0, 1:5], "3-D free surface: " * (case == "flat3D" ? "flat" :
        "oblique plane θ = $(Int(case3DTheta(case)))°, azimuth $(Int(case3DAzimuth(case)))°") *
        " -- raw displacement, one scale per line for every code; black = flat SPECFEM3D rotated (exact)"; fontsize=16)
    for (il, line) in enumerate(("deep", "shallow"))
        ks = keysOf(line)
        peak = maximum(maximum(abs, trace(reference, k, c)) for k in ks, c in (:u1, :u2, :un))
        gain = 0.45 / peak
        panels = ((:un, false), (:u1, false), (:u2, false), (:un, true), (:u1, true))
        for (ic, (comp, residual)) in enumerate(panels)
            ax = Axis(fig[il, ic]; xlabel="time (s)", yticks=(1:9, [@sprintf("%s=%g", i <= 5 ? "x'" : "y'",
                    i <= 5 ? GC3.lineOffsets[i] / 1e3 : GC3.crossOffsets[i-5] / 1e3) for i in 1:9]),
                title="$line: " * (residual ? "$(comp) − reference (× 5)" : string(comp)))
            for (ir, key) in enumerate(ks)
                ref = trace(reference, key, comp)
                residual || lines!(ax, tt, ir .+ gain .* ref; style["reference"]..., label=ir == 1 ? label["reference"] : nothing)
                for code in cs
                    v = trace(data[case][code], key, comp)
                    residual && (v = 5 .* (v .- ref))
                    lines!(ax, tt, ir .+ gain .* v; style[code]..., label=ir == 1 ? label[code] : nothing)
                end
            end
            il == 1 && ic == 1 && axislegend(ax; position=:lt, labelsize=10)
        end
        Label(fig[il, 6], @sprintf("scale:\n%.2e m\nper 0.45", peak); fontsize=11, tellheight=false)
    end
    path = joinpath(outDir, "wiggles3D_$(case).png")
    save(path, fig)
    println("figure -> ", path)
end

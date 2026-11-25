# dev/cutCells/compareVariants.jl -- tilted-plane comparison of every free-surface variant
#
# Reads the cached tiltedPlane runs (data/dev/cutCells, produced by run.jl tiltedPlane with
# FLEXOPT_CUT_SCHEME = OPT3|OPT5 and FLEXOPT_CUT_GEOMETRY = cut|zigzag) and compares, with RAW
# traces (no normalisation, no amplitude fit):
#   * each variant at θ = 10/20/30° against its own θ = 0 run (surface on grid nodes, where
#     cut and zigzag coincide) -- the exact answer is the flat solution rotated;
#   * every θ = 0 run against the cached SPECFEM2D flat run (deep line);
#   * the no-surface control (full*, OPT3) when cached.
# Figures: one per angle and receiver line, wiggles of u_n and u_t, all variants overlaid.
#
# Run (after the tiltedPlane runs):  julia --project=. dev/cutCells/compareVariants.jl

if !@isdefined(flexopt_root)
    import Pkg
    global flexopt_root = abspath(joinpath(@__DIR__, "..", ".."))
    Pkg.activate(flexopt_root)
end
!@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
using .commonBatchs
!@isdefined(specfemBenchmark) && include(joinpath(flexopt_root, "src", "specfemBenchmark.jl"))
using DrWatson, JLD2, CairoMakie, Printf, LinearAlgebra, Statistics
CairoMakie.activate!(type="png")

const CV = (dx=parse(Float64, get(ENV, "FLEXOPT_TILT_DX", "500.0")), halfWidth=60e3, duration=9.0,
    codeVersion="v1-2026-11-24", receiversX=(-30e3, -20e3, -10e3, 0.0, 10e3))
variants = [("OPT3", "cut"), ("OPT3", "zigzag"), ("OPT5", "cut"), ("OPT5", "zigzag")]
label(s, g) = g == "cut" ? "$s cut-cell" : "$s zigzag"
colors = Dict(("OPT3", "cut") => :dodgerblue, ("OPT3", "zigzag") => :darkorange,
    ("OPT5", "cut") => :crimson, ("OPT5", "zigzag") => :purple)
styles = Dict("cut" => :solid, "zigzag" => :dash)

function cachedRun(θ, scheme, geometry; freeSurface=true)
    extension = scheme == "OPT5" ? "cubic" : "linear"
    params = Dict{String,Any}("stage" => "tiltedPlane", "theta" => Float64(θ), "dx" => CV.dx,
        "halfWidth" => CV.halfWidth, "duration" => CV.duration, "extension" => extension,
        "codeVersion" => CV.codeVersion)
    freeSurface || (params["freeSurface"] = false)
    scheme == "OPT3" || (params["scheme"] = scheme)
    geometry == "cut" || (params["geometry"] = geometry)
    file = datadir("dev/cutCells", "tiltedPlane_$(commonBatchs.hash_parameters(params)).jld2")
    isfile(file) ? load(file) : nothing
end

runs = Dict{Any,Any}()
for (s, g) in variants, θ in (0, 10, 20, 30)
    r = cachedRun(θ, s, g)
    isnothing(r) || (runs[(s, g, θ)] = r)
end
for θ in (0, 10, 20, 30)
    r = cachedRun(θ, "OPT3", "cut"; freeSurface=false)
    isnothing(r) || (runs[("OPT3", "noSurface", θ)] = r)
end
println("cached runs: ", sort([string(k) for k in keys(runs)]))

function sampleAt(times, v, t)
    k = clamp(searchsortedlast(times, t), 1, length(times) - 1)
    w = clamp((t - times[k]) / (times[k+1] - times[k]), 0.0, 1.0)
    (1 - w) * v[k] + w * v[k+1]
end

println("\n== raw rel. L2 against the variant's own θ = 0 run (both components, |u| norm), window before edge returns")
@printf("  %-16s %-8s %8s %8s %8s\n", "variant", "line", "10°", "20°", "30°")
summary = Dict{Any,Float64}()
for key in vcat(variants, [("OPT3", "noSurface")]), line in ("deep", "shallow")
    s, g = key
    haskey(runs, (s, g, 0)) || continue
    ref = runs[(s, g, 0)]
    vals = String[]
    for θ in (10, 20, 30)
        haskey(runs, (s, g, θ)) || (push!(vals, "     -"); continue)
        r = runs[(s, g, θ)]
        errs = Float64[]
        for ir in eachindex(CV.receiversX)
            tr, rf = r["traces"]["$(line)_$ir"], ref["traces"]["$(line)_$ir"]
            tEnd = min(CV.duration, tr.tReturn, rf.tReturn)
            keep = ref["times"] .<= tEnd
            tt = ref["times"][keep]
            dn = [sampleAt(r["times"], tr.un, t) for t in tt] .- rf.un[keep]
            dt = [sampleAt(r["times"], tr.ut, t) for t in tt] .- rf.ut[keep]
            push!(errs, sqrt(sum(abs2, dn) + sum(abs2, dt)) / sqrt(sum(abs2, rf.un[keep]) + sum(abs2, rf.ut[keep])))
        end
        summary[(s, g, line, θ)] = mean(errs)
        push!(vals, @sprintf("%8.3f", mean(errs)))
    end
    @printf("  %-16s %-8s %s\n", g == "noSurface" ? "$s no surface" : label(s, g), line, join(vals, " "))
end

# θ = 0 against SPECFEM2D (deep line, u_z = u_n)
specfemOut = joinpath(flexopt_root, "data", "specfem2d_benchmarks", "homogeneous_flat_free_surface_dx500m_nrec23", "OUTPUT_FILES")
specfem = nothing
if CV.dx == 500.0 && isdir(specfemOut)
    f0 = 1.5 * 3464.0 / (16 * 500.0)
    stations = ["S0016", "S0022", "S0017", "S0023", "S0018"]
    files = specfemBenchmark.find_specfem2d_traces(specfemOut; component=:z, seismotype=2, network="FX", stations=stations)
    byStation = Dict(split(basename(f), '.')[2] => f for f in files)
    specfem = map(stations) do name
        raw = specfemBenchmark.read_specfem2d_trace(byStation[name]; time_shift=1.2 / f0)
        d = zeros(length(raw.values)); d[2:end] .= cumsum(((raw.values[1:end-1] .+ raw.values[2:end]) ./ 2) .* diff(raw.time))
        (time=Float64.(raw.time), values=d)
    end
    println("\n== θ = 0 (flat surface on grid nodes) vs SPECFEM2D, deep line u_z, raw rel. L2 over 0-9 s")
    for s in ("OPT3", "OPT5")
        haskey(runs, (s, "cut", 0)) || continue
        r = runs[(s, "cut", 0)]
        errs = map(eachindex(stations)) do i
            keep = specfem[i].time .<= CV.duration
            mine = [sampleAt(r["times"], r["traces"]["deep_$i"].un, t) for t in specfem[i].time[keep]]
            norm(mine .- specfem[i].values[keep]) / norm(specfem[i].values[keep])
        end
        @printf("  %-5s  %s\n", s, join([@sprintf("x'=%+3.0f km %.3f", CV.receiversX[i] / 1e3, errs[i]) for i in eachindex(errs)], "   "))
    end
end

outdir = joinpath(@__DIR__, "outputs", "tiltedPlane")
mkpath(outdir)
for θ in (10, 20, 30), line in ("deep", "shallow")
    any(haskey(runs, (s, g, θ)) for (s, g) in variants) || continue
    fig = Figure(size=(1500, 240 * length(CV.receiversX) + 90))
    Label(fig[0, 1:2], "Tilted free surface θ = $(θ)°, $(line) receivers -- raw rotated components; black: OPT5 at θ = 0 " *
        "(flat solution, the exact answer rotated)"; fontsize=14)
    reference = get(runs, ("OPT5", "cut", 0), get(runs, ("OPT3", "cut", 0), nothing))
    for (ir, xp) in enumerate(CV.receiversX), (ic, comp) in enumerate((:un, :ut))
        ax = Axis(fig[ir, ic]; ylabel="x' = $(round(Int, xp / 1e3)) km",
            title=ir == 1 ? (comp == :un ? "u_n (normal)" : "u_t (tangential)") : "",
            xlabel=ir == length(CV.receiversX) ? "time (s)" : "")
        key = "$(line)_$ir"
        isnothing(reference) || lines!(ax, reference["times"], getproperty(reference["traces"][key], comp);
            color=:black, linewidth=3, label="flat (θ = 0), OPT5")
        for (s, g) in variants
            haskey(runs, (s, g, θ)) || continue
            r = runs[(s, g, θ)]
            lines!(ax, r["times"], getproperty(r["traces"][key], comp); color=colors[(s, g)],
                linestyle=styles[g], linewidth=1.6, label=label(s, g))
        end
        ir == 1 && ic == 1 && axislegend(ax; position=:lt, nbanks=2)
    end
    path = joinpath(outdir, "variants_theta$(θ)_$(line).png")
    save(path, fig)
    println("figure -> ", path)
end

# dev/moonRSGT/run.jl -- receiver-side strain Green's tensor (RSGT) catalogue
# for a family of 1-D Moon models (.nd files), so that synthetics for any
# station, source position (down to the top of the fluid core) and
# mechanism can be made later without solving anything (lib/moonSynthetics.jl).
#
# The .nd models are CONFIDENTIAL: they are read from --nd-dir (default
# ENV MOON_ND_DIR, else /gpfs/scratch/nobuaki/moonRSGT/OUTPUT_ND) and are
# never copied into the repository; the catalogue goes to --out-dir
# (default ENV MOON_RSGT_DIR, else /gpfs/scratch/nobuaki/moonRSGT/catalogue).
#
# Steps (from the repository root):
#   julia --project=. dev/moonRSGT/run.jl prepare              # model summary, storage estimate, OPT3 recipes cached
#   julia --project=. dev/moonRSGT/run.jl check  [--model=1]   # build 1 model + compare with exact source-side synthetics
#   julia --project=. dev/moonRSGT/run.jl build --models=1:10  # (indices in the sorted .nd list, or "all")
#   julia --project=. dev/moonRSGT/run.jl build --task=7 --per-task=10   # Slurm array task 7 → models 61:70
#   julia --project=. dev/moonRSGT/run.jl info                 # which catalogue files exist
# Options: --nd-dir=DIR --out-dir=DIR --force (recompute existing files)
#   parameter overrides (defaults in lib/moonCatalogue.jl, moonDefaultParams):
#   --scheme=OPT3|DSMconv --h=5 --dz=20 --dtheta=0.5 --fmax=0.05 --tlen=3276.8 --fitTol=0.005

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))     # repo root + Pkg.activate (nothing else used)

function parseMoonArgs(args)
    pos = filter(a -> !startswith(a, "--"), args)
    step = isempty(pos) ? "prepare" : pos[1]
    step in ("prepare", "check", "build", "info") || error("unknown step '$step' (prepare, check, build, info)")
    opt = Dict{String,String}()
    for a in args
        startswith(a, "--") || continue
        kv = split(a[3:end], '='; limit=2)
        opt[kv[1]] = length(kv) == 2 ? kv[2] : "true"
    end
    known = ("nd-dir", "out-dir", "models", "task", "per-task", "model", "force",
             "scheme", "h", "dz", "dtheta", "fmax", "tlen", "fitTol")
    for k in keys(opt)
        k in known || error("unknown option --$k")
    end
    return step, opt
end

MOON_STEP, MOON_OPT = parseMoonArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS)
ndDir = get(MOON_OPT, "nd-dir", get(ENV, "MOON_ND_DIR", "/gpfs/scratch/nobuaki/moonRSGT/OUTPUT_ND"))
outDir = get(MOON_OPT, "out-dir", get(ENV, "MOON_RSGT_DIR", "/gpfs/scratch/nobuaki/moonRSGT/catalogue"))
isdir(ndDir) || error("no .nd folder $ndDir (copy OUTPUT_ND there by hand, or give --nd-dir=)")
ndFiles = sort(filter(f -> endswith(f, ".nd") && !startswith(f, "._"), readdir(ndDir)))
isempty(ndFiles) && error("no .nd files in $ndDir")
println("moonRSGT: step $MOON_STEP, $(length(ndFiles)) models in $ndDir, catalogue in $outDir")

include(joinpath(@__DIR__, "env.jl"))

params = moonDefaultParams()
for (k, key, T) in (("scheme", "scheme", String), ("h", "h", Float64), ("dz", "dz", Float64),
                    ("dtheta", "dθ", Float64), ("fmax", "fmax", Float64), ("tlen", "tlen", Float64),
                    ("fitTol", "fitTol", Float64))
    haskey(MOON_OPT, k) && (params[key] = T === String ? MOON_OPT[k] : parse(T, MOON_OPT[k]))
end
println("parameters: ", params)
catFileOf(nd) = joinpath(outDir, "rsgt_" * replace(nd, ".nd" => "") * ".jld2")
# figures only in `check`; loaded at top level (own statement: world age)
MOON_STEP == "check" && @eval using CairoMakie
MOON_STEP == "check" && CairoMakie.activate!(type="png")

"Model indices (in the sorted .nd list) selected by --models / --task --per-task."
function selectedModels()
    if haskey(MOON_OPT, "task")
        k = parse(Int, MOON_OPT["task"]); n = parse(Int, get(MOON_OPT, "per-task", "10"))
        return ((k-1)*n+1):min(k * n, length(ndFiles))
    end
    s = get(MOON_OPT, "models", "all")
    s == "all" && return 1:length(ndFiles)
    if occursin(':', s)
        a, b = parse.(Int, split(s, ':'))
        return a:min(b, length(ndFiles))
    end
    return [parse(Int, x) for x in split(s, ',')]
end

if MOON_STEP == "prepare"
    infos = [moonModelFromND(joinpath(ndDir, f); fitTol=params["fitTol"]).info for f in ndFiles]
    nz = [i.nzone for i in infos]; cmb = [i.cmbDepth for i in infos]
    moho = [length(i.layers) > 1 ? i.layers[1].dbot : NaN for i in infos]
    @printf("zones per model %d-%d, worst cubic-fit misfit %.4f; Moho %.1f-%.1f km; CMB %.1f-%.1f km\n",
            minimum(nz), maximum(nz), maximum(i.fitMisfit for i in infos),
            minimum(moho), maximum(moho), minimum(cmb), maximum(cmb))
    println("first model: ", moonModelSummary(infos[1]))
    nω = floor(Int, params["fmax"] * params["tlen"]) + 1
    nθ = length(0:params["dθ"]:180)
    nr = length(catalogueRadii(infos[1]; dz=params["dz"]).r_km)
    per = 8 * 15 * nω * nr * nθ
    @printf("catalogue: %d ω × %d r × %d θ × 15 pairs (ComplexF32) = %.2f GB per model, %.0f GB for %d models\n",
            nω, nr, nθ, per / 1e9, per * length(ndFiles) / 1e9, length(ndFiles))
    println("building the $(params["scheme"]) systems of the first model (caches the recipes in data/semiSymbolic) ...")
    mm = moonModelFromND(joinpath(ndDir, ndFiles[1]); fitTol=params["fitTol"])
    @time psv, sh = moonSystems(mm.model, params; recipes=moonRecipes, configs=shConfigs)
    println("P-SV dofs ", maximum(maximum, psv.dofmap), ", SH dofs ", maximum(maximum, sh.dofmap),
            ", zones ", length(psv.segments), ". Ready for check / build.")
elseif MOON_STEP == "check"
    k = parse(Int, get(MOON_OPT, "model", "1"))
    nd = ndFiles[k]
    @time buildMoonRSGT(joinpath(ndDir, nd), catFileOf(nd); params=params, recipes=moonRecipes,
                        configs=shConfigs, force=haskey(MOON_OPT, "force"))
    moonCheckRows = moonCheck(joinpath(ndDir, nd), catFileOf(nd), params; recipes=moonRecipes, configs=shConfigs,
                              outDir=joinpath(@__DIR__, "outputs", "check"))
elseif MOON_STEP == "build"
    sel = selectedModels()
    println("models ", first(sel), " … ", last(sel), " (", length(sel), ")")
    for (n, k) in enumerate(sel)
        nd = ndFiles[k]
        println("\n[", n, "/", length(sel), "] ", nd)
        try
            @time buildMoonRSGT(joinpath(ndDir, nd), catFileOf(nd); params=params, recipes=moonRecipes,
                                configs=shConfigs, force=haskey(MOON_OPT, "force"))
        catch err
            @error "model $nd failed" exception = (err, catch_backtrace())
        end
    end
elseif MOON_STEP == "info"
    done = [f for f in ndFiles if isfile(catFileOf(f))]
    missingModels = setdiff(ndFiles, done)
    println(length(done), " / ", length(ndFiles), " catalogue files present in ", outDir)
    isempty(missingModels) || println("missing (indices): ", join(findall(in(missingModels), ndFiles), ","))
end

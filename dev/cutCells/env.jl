# dev/cutCells/env.jl -- full flexOPT environment for the operator steps
# (freeEnd1D, tiltedPlane); the `geometry` step does not need it.

if !@isdefined(CUTCELLS_FULL_ENV)
    !@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
    @eval using .commonBatchs
    !@isdefined(backend) && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
    !@isdefined(flexOPT) && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
    @eval using .flexOPT: set_default_form!, DEFAULT_FORM, makeOPTsemiSymbolic, getModelPoints,
        numericalOperatorConstruction, prepareLinearSystem, propagateLinearSystem
    global CUTCELLS_FULL_ENV = true
end
flexOPT.set_default_form!(:weak)

# scheme of the operator steps: ENV FLEXOPT_CUT_SCHEME=OPT3 (default) | OPT5
CUT_SCHEME = get(ENV, "FLEXOPT_CUT_SCHEME", "OPT3")
CUT_SCHEME in ("OPT3", "OPT5") || error("FLEXOPT_CUT_SCHEME must be OPT3 or OPT5")
# void nodes read by cut rows: within the stencil half-width, extended with a matching degree
CUT_REACH = CUT_SCHEME == "OPT5" ? 2 : 1
# ENV FLEXOPT_CUT_EXTENSION = constant | linear | quadratic | cubic (default: linear for OPT3, cubic for OPT5)
CUT_EXTENSION = Symbol(get(ENV, "FLEXOPT_CUT_EXTENSION", CUT_SCHEME == "OPT5" ? "cubic" : "linear"))

# surface geometry: ENV FLEXOPT_CUT_GEOMETRY = cut (default, exact cut cells) | zigzag (stair-step)
CUT_GEOMETRY = get(ENV, "FLEXOPT_CUT_GEOMETRY", "cut")
CUT_GEOMETRY in ("cut", "zigzag") || error("FLEXOPT_CUT_GEOMETRY must be cut or zigzag")
CUT_VARIANT = CUT_GEOMETRY == "cut" ? CUT_SCHEME : "$(CUT_SCHEME)_zigzag"     # label in figures/files

"The OPT5 variant of an OPT3 parameter Dict, exactly as dev/freeSurfacePointSource/run.jl builds it."
function schemeParameters(opt3::Dict, scheme::AbstractString)
    scheme == "OPT3" && return opt3
    p = copy(opt3)
    p["orderBspace"] = 1
    p["pointsInSpace"] = 5
    p["fieldItpl"] = merge(opt3["fieldItpl"], (offsetSpace=2.0,))
    p["materItpl"] = merge(opt3["materItpl"], (offsetSpace=2.0,))
    return p
end

"OPT recipe cached exactly as dev/freeSurfacePointSource/run.jl does (data/semiSymbolic, version 9)."
function cachedRecipe(parameters, prefix)
    cacheParameters = Dict{String,Any}(k => v for (k, v) in parameters if k != "recipe_backend")
    cacheParameters["recipe_cache_version"] = 9
    # ENV FLEXOPT_CPU_RECIPE=1: separate cache entry built on the CPU in Float64 (otherwise the
    # existing entry is reused, whatever backend built it -- e.g. the Float32 Metal recipes of Aug. 2026)
    get(ENV, "FLEXOPT_CPU_RECIPE", "0") == "1" && (cacheParameters["recipe_backend_resolved"] = "cpu")
    myProduceOrLoad(cacheParameters, "semiSymbolic", prefix) do config
        runtime = Dict{String,Any}(config)
        pop!(runtime, "hash_id", nothing)
        pop!(runtime, "recipe_cache_version", nothing)
        runtime["recipe_backend"] = backend
        flexOPT.makeOPTsemiSymbolic(runtime)
    end
end

"Tolerance of `cutCellSelfCheck` for a recipe: Float32 recipes (GPU) agree only to Float32 precision."
function selfCheckTolerance(recipe)
    recette = recipe["recette"]
    backendType = hasproperty(recette, :recipe_backend_type) ? string(recette.recipe_backend_type) : "unknown"
    float64 = occursin("CPU", backendType)
    float64 || @warn "volume recipe built with $backendType (Float32): self-check at Float32 precision; " *
        "FLEXOPT_CPU_RECIPE=1 rebuilds it in Float64 on the CPU"
    return float64 ? 1e-8 : 2e-5
end

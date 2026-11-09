#!/usr/bin/env julia
#
# helpers/runPoissonMixedTestFunction.jl   (Claude-authored; does not modify src/)
#
# flexOPT test of mixed test functions for the Poisson equation (κ T')' = f.
#
# Derivation (helpers/optTestFunctionPrototype.jl and
# helpers/derivations/optTestFunctionMoments_sympy.py): in 1-D the truncation
# error depends on W only through its moments M_m = ∫W s^m / ∫W. The leading
# heterogeneous error term can be cancelled for every κ, T when
#   P = 5:  W = 4/15  hat₁ + 11/15  hat₂    (predicted O(h⁶) instead of O(h⁴))
#   P = 7:  W = 40/231 hat₁ + 191/231 hat₂  (predicted O(h⁸) instead of O(h⁶))
# and never for P = 3 (O(h²) for any test function once κ varies).
# hat_a is the hat on [ν-a, ν+a]; hat₂ is obtained in flexOPT through the knots
# of W: trialFunctionRefPoints = [ν-2, ν, ν+2] (e.g. [1,3,5] for P = 5).
#
# The weak form is linear in W, so the mixture is assembled from the recipes
# of its components, each normalised by ΣΓ = ∫W:
#     A = Σ_i w_i A_i / ΣΓ_i,    Γ = Σ_i w_i Γ_i / ΣΓ_i.
# In 2-D the tested W is the tensor product w(x) w(y) of the 1-D mixture
# (4 component recipes); whether that is enough in 2-D is one of the questions.
#
# Measure: boundary-free periodic residual RMS[(A T - Γ f)/Σ|Γ|], same material
# scenarios and machinery as scripts/runWavePeriodicFDvsOPT.jl.
#
# Usage (repository root):
#   FLEXOPT_TF_QUICK=1 julia --project=. helpers/runPoissonMixedTestFunction.jl
#   julia --project=. helpers/runPoissonMixedTestFunction.jl
#   FLEXOPT_TF_DIMS=2 FLEXOPT_TF_POINTS_2D=5 julia --project=. helpers/runPoissonMixedTestFunction.jl
#
# Environment: FLEXOPT_TF_DIMS (1,2)  FLEXOPT_TF_POINTS_1D (3,5,7)
#              FLEXOPT_TF_POINTS_2D (5)  FLEXOPT_TF_MAX_N_2D (32)  FLEXOPT_TF_RESUME (1)
# Outputs: scripts/tmp/poisson_mixed_test_function/

include(joinpath(@__DIR__, "..", "scripts", "runWavePeriodicFDvsOPT.jl"))

using CairoMakie
CairoMakie.activate!()

const TF_SCHEMA_VERSION = "poisson_mixed_test_function_v1"

tf_list(name, default) = String.(filter(!isempty, strip.(split(get(ENV, name, default), ","))))
tf_ints(name, default) = parse.(Int, tf_list(name, default))

const TF_QUICK = get(ENV, "FLEXOPT_TF_QUICK", "0") == "1"
const TF_RESUME = get(ENV, "FLEXOPT_TF_RESUME", "1") == "1"
const TF_DIMS = TF_QUICK ? [1] : tf_ints("FLEXOPT_TF_DIMS", "1,2")
const TF_POINTS = Dict(
    1 => TF_QUICK ? [5] : tf_ints("FLEXOPT_TF_POINTS_1D", "3,5,7"),
    2 => tf_ints("FLEXOPT_TF_POINTS_2D", "5"),
)
const TF_SIZES = Dict(
    1 => TF_QUICK ? [16, 32] : [16, 24, 32, 48, 64, 96, 128],
    2 => filter(<=(parse(Int, get(ENV, "FLEXOPT_TF_MAX_N_2D", "32"))), [8, 12, 16, 24, 32, 48]),
)
const TF_SCENARIOS = TF_QUICK ? ["homogeneous", "long_material_phase0"] : [
    "homogeneous", "same_wave_phase0", "same_wave_phase_pi2",
    "long_material_phase0", "short_material_phase0",
]
# 1-D mixture weights (hat₁, hat₂), from the moment conditions
const TF_WEIGHTS = Dict(5 => (4 / 15, 11 / 15), 7 => (40 / 231, 191 / 231))
const TF_OUTPUT_DIR = joinpath(FLEXOPT_ROOT, "scripts", "tmp", "poisson_mixed_test_function")

# ---------------------------------------------------------------- W components

"Knots of hat_a centred on the stencil centre ν (1-based node indices); `nothing` = default (all nodes)."
function hat_knots(P, a)
    a == 1 && return nothing
    c = (P + 1) ÷ 2
    first_knot = c - a * ((c - 1) ÷ a)
    return collect(first_knot:a:P)
end

"Per-axis knot specification for trialFunctionRefPoints."
function refpoints(dimension, P, spacings)
    knots = [something(hat_knots(P, a), collect(1:P)) for a in spacings]
    dimension == 1 && return all(==(1), spacings) ? nothing : knots[1]
    all(==(1), spacings) && return nothing
    return Tuple(knots)
end

"Test-function variants: name => list of (per-axis hat spacings, weight)."
function tf_variants(dimension, P)
    variants = Pair{String,Vector{Any}}[]
    push!(variants, "hat1" => Any[(ntuple(_ -> 1, dimension), 1.0)])
    P >= 5 || return variants
    push!(variants, "hat2" => Any[(ntuple(_ -> 2, dimension), 1.0)])
    if haskey(TF_WEIGHTS, P)
        w = TF_WEIGHTS[P]
        parts = Any[]
        for combo in Iterators.product(ntuple(_ -> (1, 2), dimension)...)
            weight = prod(w[a] for a in combo)
            push!(parts, (Tuple(combo), weight))
        end
        label = @sprintf("mixture %.4g hat1 + %.4g hat2", w[1], w[2])
        push!(variants, (dimension == 1 ? label : label * " (tensor)") => parts)
    end
    return variants
end

# ---------------------------------------------------------------- recipes

central_itpl(P) = (ptsSpace=1, ptsTime=1, offsetSpace=(P - 1) / 2, offsetTime=0.0,
    YorderBspace=-1, YorderBtime=-1)

function tf_component_recipe(dimension, P, spacings, dx)
    equation_id = dimension == 1 ? "1DpoissonHetero" : "2DpoissonHetero"
    params = Dict{String,Any}(
        "famousEquationType" => equation_id,
        "Δ" => dimension == 1 ? Float64(dx) : Tuple(fill(Float64(dx), dimension)),
        "orderBtime" => 0,
        "orderBspace" => 1,
        "pointsInSpace" => P,
        "pointsInTime" => 1,
        "supplementaryOrder" => 2,
        "fieldItpl" => central_itpl(P),
        "materItpl" => central_itpl(P),
        "nuGeometryMode" => :middle,
        "hierarchicalTestFunctions" => false,
        "evenOrderHalfShiftMode" => :none,
        "taylorInverseMode" => :hierarchical_constrained,
        "recipe_backend" => CPU(),
        "variationalForm" => :weak,
    )
    knots = refpoints(dimension, P, spacings)
    knots === nothing || (params["trialFunctionRefPoints"] = knots)
    return makeOPTsemiSymbolic(params)["recette"]
end

"ΣΓ of a recipe (= ∫W up to the common h scaling); placeholder force symbols are set to one."
function gamma_sum(recipe)
    compiled = flexOPT._compile_coefficient_recipes(recipe.rhs.Γjiννᶜ)
    return sum(evaluate_coefficient(compiled[j, 1, 1, 1], Dict{Any,Any}())
               for j in axes(compiled, 1))
end

function combine_recipes(recipes, weights)
    reference = first(recipes)
    for r in recipes
        size(r.lhs.Ajiννᶜ) == size(reference.lhs.Ajiννᶜ) ||
            error("component recipes have different stencils")
    end
    scales = [w / gamma_sum(r) for (r, w) in zip(recipes, weights)]
    A = sum(scale .* r.lhs.Ajiννᶜ for (r, scale) in zip(recipes, scales))
    G = sum(scale .* r.rhs.Γjiννᶜ for (r, scale) in zip(recipes, scales))
    return merge(reference, (
        lhs=merge(reference.lhs, (; Ajiννᶜ=A)),
        rhs=merge(reference.rhs, (; Γjiννᶜ=G)),
    ))
end

short_err(err) = first(replace(sprint(showerror, err), '\n' => ' '), 300)

# ---------------------------------------------------------------- run

function tf_run!(rows, done, dimension, checkpoint)
    equation_id = dimension == 1 ? "1DpoissonHetero" : "2DpoissonHetero"
    equation = only(filter(e -> e.equation == equation_id, EQUATIONS))
    scenarios = filter(s -> s.name in TF_SCENARIOS, material_scenarios(equation))
    for P in TF_POINTS[dimension], n in TF_SIZES[dimension]
        dx = 2pi / n
        components = Dict{Any,Any}()          # spacings => recipe (shared by the variants)
        for (name, parts) in tf_variants(dimension, P)
            key = (dimension, P, name, n)
            key in done && continue
            @info "mixed test function" dimension P n W=name
            status = "ok"
            numeric = nothing
            t0 = time()
            try
                recipes = map(parts) do (spacings, _)
                    get!(() -> tf_component_recipe(dimension, P, spacings, dx), components, spacings)
                end
                numeric = compile_recipe(combine_recipes(recipes, [p[2] for p in parts]))
            catch err
                err isa InterruptException && rethrow()
                status = "recipe failed: " * short_err(err)
                @warn "recipe failed" P n W=name status
            end
            seconds = time() - t0
            for scenario in scenarios
                absolute = NaN
                relative = NaN
                row_status = status
                if numeric !== nothing
                    try
                        absolute, relative, _ = residual_error(numeric, equation, scenario, n, dx, nothing)
                    catch err
                        err isa InterruptException && rethrow()
                        row_status = "residual failed: " * short_err(err)
                    end
                end
                push!(rows, (; dimension, P, W=name, scenario=scenario.name, n, dx,
                    absolute=Float64(absolute), relative=Float64(relative),
                    seconds=Float64(seconds), status=row_status))
            end
            push!(done, key)
            checkpoint()
        end
    end
end

function tf_orders(rows)
    out = NamedTuple[]
    for key in unique((r.dimension, r.P, r.W, r.scenario) for r in rows)
        rs = sort(filter(r -> (r.dimension, r.P, r.W, r.scenario) == key &&
            isfinite(r.absolute) && r.absolute > 0, rows); by=r -> r.n)
        length(rs) < 2 && continue
        orders = [log(rs[i].absolute / rs[i+1].absolute) / log(rs[i+1].n / rs[i].n)
                  for i in 1:length(rs)-1]
        push!(out, (; dimension=key[1], P=key[2], W=key[3], scenario=key[4],
            n_finest=rs[end].n, error_finest=rs[end].absolute,
            order_last=orders[end], order_median=median(orders)))
    end
    return out
end

function tf_print(summary)
    for (dimension, P) in unique((s.dimension, s.P) for s in summary)
        println("\n== $(dimension)-D, P = $P   (order at the last pair of grids; error at the finest grid)")
        sel = filter(s -> s.dimension == dimension && s.P == P, summary)
        scenarios = unique(s.scenario for s in sel)
        @printf("%-40s", "W")
        foreach(sc -> @printf("%24s", first(sc, 22)), scenarios)
        println()
        for W in unique(s.W for s in sel)
            @printf("%-40s", W)
            for sc in scenarios
                m = filter(s -> s.W == W && s.scenario == sc, sel)
                isempty(m) ? @printf("%24s", "-") :
                    @printf("%24s", @sprintf("%.2f  (%.2e)", m[1].order_last, m[1].error_finest))
            end
            println()
        end
    end
end

const TF_COLORS = Dict(3 => :gray40, 4 => :purple, 5 => :royalblue, 6 => :teal, 7 => :crimson)

function tf_plot(rows, dimension)
    sel = filter(r -> r.dimension == dimension && isfinite(r.absolute) && r.absolute > 0, rows)
    isempty(sel) && return nothing
    scenarios = unique(r.scenario for r in sel)
    fig = Figure(size=(300 * length(scenarios) + 300, 340))
    local ax_last
    for (j, scenario) in enumerate(scenarios)
        ax = Axis(fig[1, j]; title=scenario, xlabel="n", ylabel=j == 1 ? "RMS[(AT-Γf)/Σ|Γ|]" : "",
            xscale=log10, yscale=log10)
        for (P, W) in unique((r.P, r.W) for r in sel)
            rs = sort(filter(r -> r.P == P && r.W == W && r.scenario == scenario, sel); by=r -> r.n)
            isempty(rs) && continue
            style = startswith(W, "mixture") ? :solid : W == "hat2" ? :dot : :dash
            scatterlines!(ax, [r.n for r in rs], [r.absolute for r in rs];
                color=get(TF_COLORS, P, :black), linestyle=style, markersize=5,
                label="P=$P $(W)")
        end
        ax_last = ax
    end
    Legend(fig[1, length(scenarios) + 1], ax_last; labelsize=9, framevisible=false)
    file = joinpath(TF_OUTPUT_DIR, "mixed_test_function_$(dimension)D" * (TF_QUICK ? "_quick" : "") * ".png")
    save(file, fig)
    return file
end

function tf_main()
    mkpath(TF_OUTPUT_DIR)
    tag = TF_QUICK ? "_quick" : ""
    output_file = joinpath(TF_OUTPUT_DIR, "poisson_mixed_test_function$(tag).jld2")
    rows = NamedTuple[]
    if TF_RESUME && isfile(output_file)
        previous = load(output_file)
        get(previous, "schema_version", "") == TF_SCHEMA_VERSION && append!(rows, previous["rows"])
    end
    done = Set{Any}((r.dimension, r.P, r.W, r.n) for r in rows)
    checkpoint() = jldsave(output_file; schema_version=TF_SCHEMA_VERSION, rows)
    for dimension in TF_DIMS
        tf_run!(rows, done, dimension, checkpoint)
    end
    summary = tf_orders(rows)
    jldsave(output_file; schema_version=TF_SCHEMA_VERSION, rows, summary)
    for (file, data) in (("rows", rows), ("summary", summary))
        isempty(data) && continue
        open(joinpath(TF_OUTPUT_DIR, "poisson_mixed_test_function$(tag)_$(file).csv"), "w") do io
            names = collect(keys(first(data)))
            println(io, join(string.(names), ","))
            for r in data
                println(io, join((v isa AbstractString ? "\"$v\"" : string(v)
                                  for v in (getproperty(r, k) for k in names)), ","))
            end
        end
    end
    tf_print(summary)
    figures = filter(!isnothing, [tf_plot(rows, d) for d in unique(r.dimension for r in rows)])
    failures = filter(r -> r.status != "ok", rows)
    isempty(failures) || @warn "Some rows failed (status column)" count=length(failures)
    @info "mixed test-function benchmark complete" output_file figures
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && tf_main()

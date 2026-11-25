# application/Outputs_1D/app.jl -- "cooking recipes": derived quantities chosen by the user as PDE expressions
#
# The solver gives a field u; what we want is often a quantity derived from it (strain, stress,
# flux, traction...).  The user writes that quantity as a Julia function of the field and of
# material arrays, e.g.
#
#     registerOutput!("flux", (u, a) -> a * ∂x(u))                 # q = a ∂x u
#     registerOutput!("mixed", (u, a, b) -> a * ∂x(u) + b * u)
#
# and it becomes a famousEquation entry, turned into an OPT recipe by the same machinery
# (makeOPTsemiSymbolic) and the same scheme parameters as the operators.  Point values are
# obtained by projection:  q = M⁻¹ B u,  B = output recipe, M = recipe of the "value" entry c·u
# (c = 1).  M can come from the same scheme or from another one (`mass=`): the OPT schemes
# modify their mass to cancel the errors of the wave operator, which is not what an output needs.
#
# Nothing in src/ is modified: entries are new methods of flexOPT.famousEquation, as in
# dev/spherical1Delastic/lib.

using Symbolics, LinearAlgebra, Printf

# ---------------------------------------------------------------------------
# user-defined outputs -> famousEquation entries
# ---------------------------------------------------------------------------

"name => (f, expression, field, materials, nmat)"
OUTPUT_REGISTRY = @isdefined(OUTPUT_REGISTRY) ? OUTPUT_REGISTRY : Dict{String,Any}()

"""
    registerOutput!(name, f)

Declare a derived quantity q = f(u, m1, m2, ...) where u is the field and m1, m2, ... material
arrays (functions of x), written with ∂x, ∂x², +, *, ... (Symbolics).  The famousEquation entry
"out_<name>" is created; its recipe cache key contains the expression, so changing f under the
same name never reuses an old recipe.  Call it as its own top-level statement (REPL line),
before using the output: the new method must be visible to the recipe builder (world age).
"""
function registerOutput!(name::AbstractString, f::Function)
    nmat = first(methods(f)).nargs - 2
    nmat <= 6 || error("at most 6 material arrays")
    @variables u(x) m1(x) m2(x) m3(x) m4(x) m5(x) m6(x)
    mats = [m1, m2, m3, m4, m5, m6][1:nmat]
    expr = Base.invokelatest(f, u, mats...)
    OUTPUT_REGISTRY[name] = (f=f, expr=expr, field=u, materials=mats, nmat=nmat)
    entry = Symbol("eq_out_" * name)
    @eval function flexOPT.famousEquation(::Val{$(QuoteNode(entry))})
        r = OUTPUT_REGISTRY[$name]
        @variables g(x)
        vars = length(r.materials) == 1 ? r.materials[1] : Tuple(r.materials)
        coordinates = (x,)
        ∂, ∂² = usefulPartials(coordinates)
        return r.expr, r.field, vars, g, g, g, coordinates, ∂, ∂²
    end
    println("output \"$name\": q = ", expr)
    return OUTPUT_REGISTRY[name]
end

# the projection "mass": c·u with c = 1
registerOutput!("value", (u, c) -> c * u)

# ---------------------------------------------------------------------------
# schemes (same parameters as application/Poisson_1D and dev/spherical1Delastic)
# ---------------------------------------------------------------------------

_itpl(off) = (ptsSpace=1, ptsTime=1, offsetSpace=off, offsetTime=1, YorderBspace=-1, YorderBtime=1)
OUTPUTS1D_SCHEMES = Dict{String,NamedTuple}(
    "convFD3" => (orderBtime=1, orderBspace=-1, pointsInSpace=3, pointsInTime=1, supplementaryOrder=0, fieldItpl=_itpl(1.0), materItpl=_itpl(1.0)),
    "convFD5" => (orderBtime=1, orderBspace=-1, pointsInSpace=5, pointsInTime=1, supplementaryOrder=0, fieldItpl=_itpl(2.0), materItpl=_itpl(2.0)),
    "OPT3"    => (orderBtime=1, orderBspace=1, pointsInSpace=3, pointsInTime=1, supplementaryOrder=2, fieldItpl=_itpl(1.0), materItpl=_itpl(1.0)),
    "OPT5"    => (orderBtime=1, orderBspace=1, pointsInSpace=5, pointsInTime=1, supplementaryOrder=2, fieldItpl=_itpl(2.0), materItpl=_itpl(2.0)),
    # plain Petrov–Galerkin (no OPT error-cancelling correction): candidate masses for outputs (experimental)
    "Gal3"    => (orderBtime=1, orderBspace=1, pointsInSpace=3, pointsInTime=1, supplementaryOrder=0, fieldItpl=_itpl(1.0), materItpl=_itpl(1.0)),
    "Gal5"    => (orderBtime=1, orderBspace=1, pointsInSpace=5, pointsInTime=1, supplementaryOrder=0, fieldItpl=_itpl(2.0), materItpl=_itpl(2.0)),
)
OUTPUTS1D_VERSION = "v2-2026-11-25"

"OPT recipe of output `name` for `scheme` at spacing Δ (cached in data/semiSymbolic)."
function outputRecipe(name::AbstractString, scheme::AbstractString, Δ::Real)
    haskey(OUTPUT_REGISTRY, name) || error("unknown output $name: registerOutput! it first")
    s = OUTPUTS1D_SCHEMES[scheme]
    p = Dict{String,Any}(
        "famousEquationType" => "out_" * name, "Δ" => Float64(Δ),
        "orderBtime" => s.orderBtime, "orderBspace" => s.orderBspace,
        "pointsInSpace" => s.pointsInSpace, "pointsInTime" => s.pointsInTime,
        "supplementaryOrder" => s.supplementaryOrder, "fieldItpl" => s.fieldItpl, "materItpl" => s.materItpl,
        "variationalForm" => s.orderBspace == -1 ? :strong : :weak,
        "outputExpression" => string(OUTPUT_REGISTRY[name].expr), "outputsVersion" => OUTPUTS1D_VERSION,
        "recipe_backend" => :cpu)
    return myProduceOrLoad(makeOPTsemiSymbolic, p, "semiSymbolic")
end

"Dense matrix of a recipe on the grid for the material arrays `mats` (path of dev/spherical1Delastic flexOPTSHMatrix)."
function outputMatrix(optRec, mats::Vector, Δ::Real; modelName="outputs1D")
    Nx = length(mats[1])
    modelFam = (models=[collect(Float64, m) for m in mats], modelName=modelName, modelPoints=(Nx), Δ=(Float64(Δ)))
    params = Dict{String,Any}("optRec" => optRec, "modelFam" => modelFam,
                              "absorbingBoundaries" => nothing, "maskedRegionInSpace" => nothing,
                              "representation" => "matrixfree")
    prepared = prepareLinearSystem(numericalOperatorConstruction(params)["numOperators"])
    knownField = copy(prepared.known_lhs_template)
    knownForce = copy(prepared.known_rhs_template); knownForce .= 0.0
    A, _ = evaluateLinearSystem(prepared, knownField, knownForce; sparse_output=false)
    return Matrix(A)
end

"""
    cookOutput(name, uvals, mats, Δ; scheme="OPT3", mass=scheme, boundary=nothing) -> (lumped, consistent)

Point values of output `name` from the nodal field `uvals` and material arrays `mats` (in the
order of the function given to registerOutput!), on a uniform grid of spacing Δ.
`boundary` = values imposed on the first/last `pointsInSpace` nodes (the recipe rows there are
not output rows yet: traction / one-sided rows are the next step); default: the lumped values.
"""
function cookOutput(name::AbstractString, uvals::AbstractVector, mats::Vector, Δ::Real;
                    scheme::AbstractString="OPT3", mass::AbstractString=scheme, boundary=nothing)
    Nx = length(uvals)
    tag = "$(name)_$(scheme)_$(Nx)_$(round(Δ, sigdigits=6))"
    B = outputMatrix(outputRecipe(name, scheme, Δ), mats, Δ; modelName="out_B_" * tag)
    M = outputMatrix(outputRecipe("value", mass, Δ), [ones(Nx)], Δ; modelName="out_M_$(mass)_" * tag)
    Bu = B * uvals
    lumped = Bu ./ (M * ones(Nx))
    p = max(OUTPUTS1D_SCHEMES[scheme].pointsInSpace, OUTPUTS1D_SCHEMES[mass].pointsInSpace)
    ends = vcat(1:p, Nx-p+1:Nx)
    bvals = isnothing(boundary) ? lumped[ends] : boundary[ends]
    for (k, i) in enumerate(ends)
        M[i, :] .= 0.0; M[i, i] = 1.0; Bu[i] = bvals[k]
    end
    return lumped, M \ Bu
end

# ---------------------------------------------------------------------------
# test against an exact answer
# ---------------------------------------------------------------------------

"Exact q(x) at the nodes for u = uex(x), materials matex[i](x) (Julia functions of a symbolic x)."
function exactOutput(name::AbstractString, uex::Function, matex::Vector, X::AbstractVector)
    f = OUTPUT_REGISTRY[name].f
    e = expand_derivatives(Base.invokelatest(f, uex(x), [m(x) for m in matex]...))
    return [Float64(Symbolics.value(substitute(e, Dict(x => xi)))) for xi in X]
end

OUTPUTS1D_TEST_U = ξ -> cos(ξ)
OUTPUTS1D_TEST_MATS = [ξ -> 2 + sin(ξ / 3), ξ -> 1 + cos(ξ / 5) / 2, ξ -> 1.5 + 0 * ξ]

"""
    outputs1DCompare(name="flux"; G=8, L=60, pairs=[...]) -> NamedTuple

u = cos x on [0, L] with G points per wavelength; materials 2 + sin(x/3), 1 + cos(x/5)/2, 1.5
(as many as the output needs).  `pairs` = (B scheme, M scheme) combinations; each gives a
"lumped" and a "consistent" curve.  Errors on the interior (12 nodes from each end excluded);
a pair whose recipe fails is reported and skipped.
"""
function outputs1DCompare(name::AbstractString="flux"; G::Real=8.0, L::Real=60.0,
                          pairs=[("convFD3", "convFD3"), ("convFD5", "convFD5"), ("OPT3", "OPT3"), ("OPT5", "OPT5"),
                                 ("OPT3", "Gal3"), ("OPT5", "Gal5"), ("Gal3", "Gal3")],
                          verbose::Bool=true)
    r = OUTPUT_REGISTRY[name]
    h = 2π / G; Nx = round(Int, L / h) + 1; h = L / (Nx - 1)
    X = collect(range(0.0, L; length=Nx))
    matex = OUTPUTS1D_TEST_MATS[1:r.nmat]
    u = OUTPUTS1D_TEST_U.(X)
    mats = [m.(X) for m in matex]
    q = exactOutput(name, OUTPUTS1D_TEST_U, matex, X)
    interior = 13:(Nx - 12)
    out = Dict{String,Vector{Float64}}()
    for (bs, ms) in pairs
        lab = bs == ms ? bs : "$(bs) + $(ms) mass"
        try
            l, c = cookOutput(name, u, mats, h; scheme=bs, mass=ms, boundary=q)
            out["$(lab) lumped"] = l; out["$(lab) consistent"] = c
        catch err
            @warn "pair ($bs, $ms) failed" exception = err
        end
    end
    err = Dict(m => norm(v[interior] .- q[interior]) / norm(q[interior]) for (m, v) in out)
    if verbose
        @printf("output %s: q = %s;  G = %.2f points/wavelength (kh = %.3f), Nx = %d\n", name, string(r.expr), G, h, Nx)
        for m in sort(collect(keys(err)); by=k -> err[k])
            @printf("   %-28s %.3e\n", m, err[m])
        end
    end
    return (; name, X, u, q, out, err, interior, G, h, Nx)
end

"Error versus points per wavelength."
function outputs1DConvergence(name::AbstractString="flux"; Gs=[4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0], kwargs...)
    rows = [outputs1DCompare(name; G=G, kwargs...) for G in Gs]
    methods = sort(collect(intersect([Set(keys(r.err)) for r in rows]...)))
    return (; name, Gs=collect(Float64, Gs), methods, err=Dict(m => [r.err[m] for r in rows] for m in methods), rows)
end

# the two outputs of the first test
registerOutput!("flux", (u, a) -> a * ∂x(u))               # 1-D analogue of σ = μ ∂x u
registerOutput!("mixed", (u, a, b) -> a * ∂x(u) + b * u)    # several terms, mixed orders

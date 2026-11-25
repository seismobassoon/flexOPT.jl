# application/Outputs_1D/app.jl -- "cooking recipes": derived quantities declared like famous equations
#
# A field u is not always what we want: strain, stress, flux, traction... are
# derived from it by a PDE expression.  Here the derived quantity is declared
# exactly like an equation (a famousEquation entry) and turned into an OPT
# recipe with the SAME machinery and scheme parameters as the operator:
#
#   output entry   "1DoutFlux"   q = a(x) ∂x u        (1-D analogue of σ = μ ∂x u)
#   mass entry     "1DoutValue"  m = c(x) u            (c = 1: the test-function "mass")
#
# For a weak-form scheme (OPT3, OPT5) the recipe rows give ∫ φ_ν q dx, not
# q(x_ν).  Two ways to get point values out of them are compared:
#   lumped      q_ν = (B u)_ν / (M 1)_ν
#   consistent  q   = M⁻¹ B u          (Petrov–Galerkin projection, both B and M from recipes)
# against the conventional Fornberg derivatives (FD3, FD5) that
# dev/spherical1DelasticKernel currently uses for strains.
#
# Test: u = cos(kx), a(x) = 2 + sin(x/3) on [0, L]; exact q = -a k sin(kx).
# Nothing in src/ is modified: the two entries are new methods of
# flexOPT.famousEquation, as in dev/spherical1Delastic/lib.

using Symbolics, LinearAlgebra, Printf

# ---------------------------------------------------------------------------
# the two "cooking" entries (coordinate x from commonBatchs)
# ---------------------------------------------------------------------------

function flexOPT.famousEquation(::Val{:eq_1DoutFlux})
    @variables u(x) a(x) g(x)
    exprs = a * ∂x(u)
    fields = u
    vars = a
    extexprs = g; extfields = g; extvars = g
    coordinates = (x,)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function flexOPT.famousEquation(::Val{:eq_1DoutValue})
    @variables u(x) c(x) g(x)
    exprs = c * u
    fields = u
    vars = c
    extexprs = g; extfields = g; extvars = g
    coordinates = (x,)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

# bump when an entry above changes (part of the recipe cache key)
const OUTPUTS1D_ENTRY_VERSION = "v1-2026-11-25"

# ---------------------------------------------------------------------------
# schemes (same parameters as application/Poisson_1D and dev/spherical1Delastic)
# ---------------------------------------------------------------------------

OUTPUTS1D_SCHEMES = Dict{String,NamedTuple}(
    "convFD3" => (orderBtime=1, orderBspace=-1, pointsInSpace=3, pointsInTime=1, supplementaryOrder=0,
        fieldItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=1),
        materItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=1)),
    "OPT3" => (orderBtime=1, orderBspace=1, pointsInSpace=3, pointsInTime=1, supplementaryOrder=2,
        fieldItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=1),
        materItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=1)),
    "OPT5" => (orderBtime=1, orderBspace=1, pointsInSpace=5, pointsInTime=1, supplementaryOrder=2,
        fieldItpl=(ptsSpace=1, ptsTime=1, offsetSpace=2.0, offsetTime=1, YorderBspace=-1, YorderBtime=1),
        materItpl=(ptsSpace=1, ptsTime=1, offsetSpace=2.0, offsetTime=1, YorderBspace=-1, YorderBtime=1)),
)

"OPT recipe of `entry` for `scheme`, built once at Δ = 1 (scaled coordinate s = x/h) and cached in data/semiSymbolic."
function outputs1DRecipe(entry::AbstractString, scheme::AbstractString)
    haskey(OUTPUTS1D_SCHEMES, scheme) || error("unknown scheme $scheme; available: $(join(sort(collect(keys(OUTPUTS1D_SCHEMES))), ", "))")
    s = OUTPUTS1D_SCHEMES[scheme]
    p = Dict{String,Any}(
        "famousEquationType" => entry, "Δ" => 1.0,
        "orderBtime" => s.orderBtime, "orderBspace" => s.orderBspace,
        "pointsInSpace" => s.pointsInSpace, "pointsInTime" => s.pointsInTime,
        "supplementaryOrder" => s.supplementaryOrder,
        "fieldItpl" => s.fieldItpl, "materItpl" => s.materItpl,
        "variationalForm" => s.orderBspace == -1 ? :strong : :weak,
        "outputsEntryVersion" => OUTPUTS1D_ENTRY_VERSION,
        "recipe_backend" => :cpu)
    return myProduceOrLoad(makeOPTsemiSymbolic, p, "semiSymbolic")
end

"Dense matrix of a recipe on Nx nodes for the material array `mat` (same path as dev/spherical1Delastic flexOPTSHMatrix)."
function outputs1DMatrix(optRec, mat::AbstractVector; modelName="outputs1D")
    Nx = length(mat)
    modelFam = (models=[collect(Float64, mat)], modelName=modelName, modelPoints=(Nx), Δ=(1.0))
    params = Dict{String,Any}("optRec" => optRec, "modelFam" => modelFam,
                              "absorbingBoundaries" => nothing, "maskedRegionInSpace" => nothing,
                              "representation" => "matrixfree")
    numOps = numericalOperatorConstruction(params)["numOperators"]
    prepared = prepareLinearSystem(numOps)
    knownField = copy(prepared.known_lhs_template)
    knownForce = copy(prepared.known_rhs_template); knownForce .= 0.0
    A, _ = evaluateLinearSystem(prepared, knownField, knownForce; sparse_output=false)
    return Matrix(A)
end

# ---------------------------------------------------------------------------
# conventional derivatives (what dev/spherical1DelasticKernel uses today)
# ---------------------------------------------------------------------------

"Fornberg (1988) weights for the m-th derivative at x0 on nodes x."
function outputs1DFornberg(x0::Real, x::AbstractVector, m::Integer=1)
    n = length(x); c = zeros(n, m + 1)
    c1 = 1.0; c4 = x[1] - x0; c[1, 1] = 1.0
    for i in 2:n
        mn = min(i - 1, m); c2 = 1.0; c5 = c4; c4 = x[i] - x0
        for j in 1:i-1
            c3 = x[i] - x[j]; c2 *= c3
            if j == i - 1
                for k in mn:-1:1
                    c[i, k+1] = c1 * (k * c[i-1, k] - c5 * c[i-1, k+1]) / c2
                end
                c[i, 1] = -c1 * c5 * c[i-1, 1] / c2
            end
            for k in mn:-1:1
                c[j, k+1] = (c4 * c[j, k+1] - k * c[j, k]) / c3
            end
            c[j, 1] = c4 * c[j, 1] / c3
        end
        c1 = c2
    end
    return c[:, m+1]
end

"Derivative matrix in s (unit spacing) with `npts`-point stencils, centred where possible, one-sided at the ends."
function outputs1DFDMatrix(Nx::Integer, npts::Integer)
    D = zeros(Nx, Nx); s = collect(0.0:Nx-1); hw = npts ÷ 2
    for i in 1:Nx
        lo = clamp(i - hw, 1, Nx - npts + 1); idx = lo:lo+npts-1
        D[i, idx] = outputs1DFornberg(s[i], s[idx], 1)
    end
    return D
end

# ---------------------------------------------------------------------------
# one grid: every method on u = cos(kx)
# ---------------------------------------------------------------------------

outputs1DMaterial(x) = 2.0 + sin(x / 3)

"""
    outputs1DCompare(; G=8.0, L=60.0, k=1.0, schemes=["convFD3","OPT3","OPT5"]) -> NamedTuple

Grid with G points per wavelength (h = 2π/(kG)), all methods, exact q.
`interior` = nodes more than 12 nodes from the ends (boundary rows of the
recipes are replaced by one-sided FD5 rows and excluded from the error).
"""
function outputs1DCompare(; G::Real=8.0, L::Real=60.0, k::Real=1.0,
                          schemes=["convFD3", "OPT3", "OPT5"], verbose::Bool=true)
    h = 2π / (k * G)
    Nx = round(Int, L / h) + 1
    h = L / (Nx - 1)
    X = collect(0.0:h:L)[1:Nx]
    a = outputs1DMaterial.(X)
    u = cos.(k .* X)
    q = -a .* k .* sin.(k .* X)
    nb = 12
    interior = (nb + 1):(Nx - nb)
    D5 = outputs1DFDMatrix(Nx, 5)
    out = Dict{String,Vector{Float64}}()
    out["FD3 (Fornberg)"] = a .* (outputs1DFDMatrix(Nx, 3) * u) ./ h
    out["FD5 (Fornberg)"] = a .* (D5 * u) ./ h
    for sc in schemes
        B = outputs1DMatrix(outputs1DRecipe("1DoutFlux", sc), a; modelName="outputs1D_flux_$(sc)_$(Nx)")
        M = outputs1DMatrix(outputs1DRecipe("1DoutValue", sc), ones(Nx); modelName="outputs1D_mass_$(sc)_$(Nx)")
        p = OUTPUTS1D_SCHEMES[sc].pointsInSpace
        # recipe rows at the ends: replaced by q_i = a_i (D5 u)_i / h (identity mass row)
        for i in vcat(1:p, Nx-p+1:Nx)
            B[i, :] .= a[i] .* D5[i, :]; M[i, :] .= 0.0; M[i, i] = 1.0
        end
        Bu = (B * u) ./ h
        out["$(sc) lumped"] = Bu ./ (M * ones(Nx))
        out["$(sc) consistent"] = M \ Bu
    end
    err = Dict(m => norm(v[interior] .- q[interior]) / norm(q[interior]) for (m, v) in out)
    if verbose
        @printf("G = %5.2f points/wavelength (kh = %.3f), Nx = %d\n", G, k * h, Nx)
        for m in sort(collect(keys(err)))
            @printf("   %-22s relative L2 error %.3e\n", m, err[m])
        end
    end
    return (; X, a, u, q, out, err, interior, G, h, Nx)
end

"Error versus points per wavelength for every method."
function outputs1DConvergence(; Gs=[4.0, 5.0, 6.0, 8.0, 10.0, 14.0, 20.0, 28.0, 40.0], kwargs...)
    rows = [outputs1DCompare(; G=G, kwargs...) for G in Gs]
    methods = sort(collect(keys(rows[1].err)))
    return (; Gs=collect(Gs), methods, err=Dict(m => [r.err[m] for r in rows] for m in methods), rows)
end

# src/cutCells/cutCellOperators.jl -- stage 2: OPT rows of cut nodes
#
# For a node whose test-function support is cut by the interface, the OPT
# row is the same weak-form integral as in the volume recipe, taken over the
# material part of the support only. Everything except the integration
# domain is the volume recipe's own:
#
#   * the Taylor inversions Cˡη, Cˡηᶜ (same `taylorCoefInv` cache),
#   * the α tables (bigαFinder) and the natural-weak-form groups,
#   * the 1-D integrands W·Yμ·Yμᶜ·Kμ·Kμᶜ of every coordinate (same
#     `WYYKKIntegralSymbolic` cache; the antiderivatives are differentiated
#     back and expanded, cell by cell, in s = 2y - 1).
#
# What changes is the final product over coordinates (`windowContraction!`
# multiplies 1-D integrals): here the spatial part is
#
#     Σ_cells Σ_a  mS_cell[a] · Π_d c_d,cell(a_d)
#
# with mS the cut-cell moments of CutCellGeometry in the s-basis. Time keeps
# its ordinary 1-D integral. Material values are evaluated numerically, so
# no symbolic recipe is built per node.
#
# Copies of the logic of `_constructAmatrix_single` (makeOPTrecette.jl) and
# of the row/column layout of `numericalOperatorConstruction`
# (makeCostFunctions.jl) and `_split_operator_columns` (linearSolution.jl);
# none of those files is modified. `cutCellSelfCheck` guards against the two
# drifting apart: a node with an uncut support must reproduce the volume row.

using Symbolics

"""
    CutCellBlock

One l.h.s. weak-form group (or the r.h.s.) of an OPT recipe, in numerical form.
"""
struct CutCellBlock
    alphaNodes::Vector{Any}                 # α.node, in the order of the α tables
    alphaExpr::Vector{Int}
    alphaField::Vector{Int}
    vars::Vector{Any}                       # material variables of the block
    table::Array{Int32,3}                   # (2 + 2nCoord, nLoop, nα), as tableForLoop
    tableμ::Matrix{Int32}                   # (nCoord, Pμ)
    tableμᶜ::Matrix{Int32}                  # (nCoord, Pμᶜ)
    C::Array{Float64,3}                     # (P, L, Pμ)
    Cc::Array{Float64,3}                    # (P, L, Pμᶜ)
    spacePieces::Vector{Array{Float64,6}}   # per space axis: (l_n+1, lᶜ_nᶜ+1, Pμ, Pμᶜ, cell, power+1), × Δ
    cellOffsets::Vector{Vector{Int}}        # per space axis: lower node of each cell relative to ν
    timeIntegral::Array{Float64,4}          # (l_n+1, lᶜ_nᶜ+1, Pμ, Pμᶜ)
    NExpr::Int
    NField::Int
end

"""
    CutCellScheme(params)

Numerical ingredients of the OPT recipe `params` (the same Dict given to
`makeOPTsemiSymbolic`; weak form, one ν geometry, no hierarchical test
functions, no clipped bounds) for cut-cell rows.
"""
struct CutCellScheme
    params::Dict
    lhs::Vector{CutCellBlock}
    rhs::CutCellBlock
    pointsShape::NTuple                     # space-time stencil shape
    pointsIndices::Array                    # SVector of each stencil point (space..., time)
    ν::Vector{Int}                          # stencil centre (space..., time)
    nSpace::Int
    maxDegree::Int                          # s-degree needed per axis
end

function CutCellScheme(params::Dict)
    @unpack famousEquationType, Δ, orderBtime, orderBspace, pointsInSpace, pointsInTime,
        supplementaryOrder, fieldItpl, materItpl = params
    for (key, default) in (("trialFunctionRefPoints", nothing), ("testIntegrationBounds", nothing),
                           ("hierarchicalTestFunctions", false), ("evenOrderHalfShiftMode", :none),
                           ("nuGeometryMode", :middle))
        value = _opt_paramget(params, Symbol(key), default)
        (isnothing(default) ? isnothing(value) : Symbol(string(value)) == Symbol(string(default))) ||
            throw(ArgumentError("cut-cell rows need the plain volume recipe: $key = $value is not supported"))
    end
    nuCentre = _opt_paramget(params, :nuCentre, nothing)
    exactTaylorTotalDegree = _opt_paramget(params, :exactTaylorTotalDegree, nothing)
    taylorInverseMode = Symbol(_opt_paramget(params, :taylorInverseMode,
        _opt_paramget(params, :taylor_inverse_mode, :hierarchical_constrained)))
    variationalForm = _normalise_variational_form(_opt_paramget(params, :variationalForm, DEFAULT_FORM[]))
    variationalForm === :natural_weak ||
        throw(ArgumentError("cut-cell rows are defined for the weak form only"))
    Δnum = SVector(Δ)
    trial = (orderBtime=orderBtime, orderBspace=orderBspace, pointsInSpace=pointsInSpace,
        pointsInTime=pointsInTime, nuGeometryMode=:middle, nuCentre=nuCentre,
        exactTaylorTotalDegree=exactTaylorTotalDegree)
    TOμ = TaylorOptions(fieldItpl, supplementaryOrder)
    TOμᶜ = TaylorOptions(materItpl, supplementaryOrder)
    eqL, eqR = famousEquations(famousEquationType)
    groups = weakTermGroups(naturalWeakForm(eqL))
    lhs = CutCellBlock[]
    geometryRef = nothing
    for orders in sort!(collect(keys(groups)); by=identity)
        eq = merge(eqL, (exprs=Tuple(groups[orders]),))
        block, geometryHere = _cutCellBlock(eq, trial, TOμ, TOμᶜ, Δnum, taylorInverseMode, collect(orders))
        isnothing(geometryRef) && (geometryRef = geometryHere)
        geometryHere.pointsShape == geometryRef.pointsShape || error("weak-form groups use different stencils")
        push!(lhs, block)
    end
    rhs, geometryR = _cutCellBlock(eqR, trial, TOμ, TOμᶜ, Δnum, taylorInverseMode,
        zeros(Int, length(eqR.coordinates)))
    geometryR.pointsShape == geometryRef.pointsShape || error("l.h.s. and r.h.s. stencils differ")
    nSpace = length(eqL.coordinates) - 1
    maxDegree = maximum(size(p, 6) - 1 for b in vcat(lhs, rhs) for p in b.spacePieces)
    CutCellScheme(params, lhs, rhs, geometryRef.pointsShape, geometryRef.pointsIndices,
        geometryRef.ν, nSpace, maxDegree)
end

function _cutCellBlock(eq, trial, TOμ, TOμᶜ, Δnum, taylor_inverse_mode, testDerivativeOrders)
    numbers = numbersOfTheExpression(eq, trial, TOμ, TOμᶜ)
    numbers.timeMarching || throw(ArgumentError("cut-cell rows are implemented for time-marching recipes"))
    _, ordersμ0, configsTaylorμ, ordersμᶜ0, configsTaylorμᶜ = investigateDependencies(eq, numbers, trial, TOμ, TOμᶜ)
    bigα, _, _ = bigαFinder(eq, numbers, ordersμ0)
    level = maximum(ordersμ0.orderBspline)
    ordersμ = _orders_at_test_level(ordersμ0, level)
    ordersμᶜ = _orders_at_test_level(ordersμᶜ0, level)
    configsTaylorμ.numberGeometries == 1 || throw(ArgumentError("one ν geometry expected"))

    # -- as _constructAmatrix_single (makeOPTrecette.jl), iConfigGeometry = 1
    @unpack nCoordinates, NtypeofExpr, NtypeofFields = numbers
    vars = eq.vars isa Tuple ? collect(eq.vars) : [eq.vars]
    multiOrdersIndices = configsTaylorμ.multiOrdersIndices
    exactTaylorTotalDegree = configsTaylorμ.exactTaylorTotalDegree
    pointsIndices = configsTaylorμ.availablePointsConfigurations[1]
    middleLinearν = configsTaylorμ.centrePointConfigurations[1]
    μPoints = configsTaylorμ.availableμPoints[1]
    μᶜPoints = configsTaylorμᶜ.availableμPoints[1]
    μaxes = configsTaylorμ.availableμaxes[1]
    μᶜaxes = configsTaylorμᶜ.availableμaxes[1]
    pointν = pointsIndices[middleLinearν]
    orderBspline = ordersμ.orderBspline
    YorderBsplineμ = ordersμ.YorderBspline
    YorderBsplineμᶜ = ordersμᶜ.YorderBspline
    L_MINUS_N = multiOrdersIndices
    L_MINUS_N = L_MINUS_N .- L_MINUS_N[1]

    timeIntegral = Array{Float64,4}(undef, 0, 0, 0, 0)
    spacePieces = Array{Float64,6}[]
    cellOffsets = Vector{Int}[]
    for iCoord in 1:nCoordinates
        maxNodes = pointsIndices[end][iCoord]
        ν = pointν[iCoord]
        lᶜ_nᶜ_max = L_MINUS_N[end][iCoord]
        l_n_max = L_MINUS_N[end][iCoord]
        νRefPoints = collect(1:maxNodes)
        params1D = (orderBspline1D=orderBspline[iCoord], YorderBspline1Dμᶜ=YorderBsplineμᶜ[iCoord],
            YorderBspline1Dμ=YorderBsplineμ[iCoord], μᶜs=μᶜaxes[iCoord], μs=μaxes[iCoord],
            maxNode=pointsIndices[end][iCoord], ν=ν, νRefPoints=νRefPoints, WνOffset=0.0,
            WDerivativeOrder=testDerivativeOrders[iCoord], integrationBounds=nothing,
            lᶜ_nᶜ_max=lᶜ_nᶜ_max, l_n_max=l_n_max, Δ=Δnum[iCoord], ImakeReport=true)
        if iCoord == nCoordinates
            timeIntegral = WYYKKIntegralNumerical(params1D)[:, :, :, :, 1]
        else
            pieces, offsets = _cutCellPieces1D(params1D)
            push!(spacePieces, pieces)
            push!(cellOffsets, offsets)
        end
    end

    coefInversionDict = Dict{String,Any}(@strdict multiOrdersIndices pointsIndices μpointsIndices=μPoints Δ=Δnum)
    coefInversionDict["pinv_version"] = "$(taylor_inverse_mode)_hierarchical_v1"
    coefInversionDict["taylor_inverse_mode"] = taylor_inverse_mode
    coefInversionDict["exact_taylor_total_degree"] = exactTaylorTotalDegree
    C = myProduceOrLoad(TaylorCoefInversion, coefInversionDict, "taylorCoefInv")["CˡηGlobal"]
    coefInversionDict = Dict{String,Any}(@strdict multiOrdersIndices pointsIndices μpointsIndices=μᶜPoints Δ=Δnum)
    coefInversionDict["pinv_version"] = "$(taylor_inverse_mode)_hierarchical_v1"
    coefInversionDict["taylor_inverse_mode"] = taylor_inverse_mode
    coefInversionDict["exact_taylor_total_degree"] = exactTaylorTotalDegree
    Cc = myProduceOrLoad(TaylorCoefInversion, coefInversionDict, "taylorCoefInv")["CˡηGlobal"]

    nLs = length(multiOrdersIndices)
    nTotalSmallα = sum(length(bigα[iExpr, iField]) for iField in 1:NtypeofFields, iExpr in 1:NtypeofExpr)
    LI_L_MINUS_N_plus_1 = LinearIndices(L_MINUS_N .+ vec2car(ones(Int, nCoordinates)))
    tableForLoop = zeros(Int32, 2 + nCoordinates * 2, nLs * nLs, nTotalSmallα)
    alphaNodes, alphaExpr, alphaField = Any[], Int[], Int[]
    indexLinearα = 1
    for iExpr in 1:NtypeofExpr, iField in 1:NtypeofFields
        for eachα in bigα[iExpr, iField]
            nᶜ = eachα.nᶜ - vec2car(ones(Int, nCoordinates))
            n = eachα.n - vec2car(ones(Int, nCoordinates))
            Lᶜ_avail = (nᶜ .+ L_MINUS_N) ∩ L_MINUS_N
            L_avail = (n .+ L_MINUS_N) ∩ L_MINUS_N
            iL = 1
            for l in L_avail, lᶜ in Lᶜ_avail
                tableForLoop[1, iL, indexLinearα] = LI_L_MINUS_N_plus_1[lᶜ+vec2car(ones(Int, nCoordinates))]
                tableForLoop[2, iL, indexLinearα] = LI_L_MINUS_N_plus_1[l+vec2car(ones(Int, nCoordinates))]
                tmplᶜ_nᶜ = lᶜ - nᶜ + vec2car(ones(Int, nCoordinates))
                tmpl_n = l - n + vec2car(ones(Int, nCoordinates))
                for iCoord in 1:nCoordinates
                    tableForLoop[2+iCoord, iL, indexLinearα] = tmplᶜ_nᶜ[iCoord]
                    tableForLoop[2+iCoord+nCoordinates, iL, indexLinearα] = tmpl_n[iCoord]
                end
                iL += 1
            end
            push!(alphaNodes, eachα.node)
            push!(alphaExpr, iExpr)
            push!(alphaField, iField)
            indexLinearα += 1
        end
    end
    tableμ = Array{Int32,2}(undef, nCoordinates, length(μPoints))
    for iμ in CartesianIndices(μPoints), iCoord in 1:nCoordinates
        tableμ[iCoord, LinearIndices(μPoints)[iμ]] = iμ[iCoord]
    end
    tableμᶜ = Array{Int32,2}(undef, nCoordinates, length(μᶜPoints))
    for iμᶜ in CartesianIndices(μᶜPoints), iCoord in 1:nCoordinates
        tableμᶜ[iCoord, LinearIndices(μᶜPoints)[iμᶜ]] = iμᶜ[iCoord]
    end

    block = CutCellBlock(alphaNodes, alphaExpr, alphaField, vars, tableForLoop, tableμ, tableμᶜ,
        Float64.(C), Float64.(Cc), spacePieces, cellOffsets, timeIntegral, NtypeofExpr, NtypeofFields)
    geometry = (pointsShape=size(pointsIndices), pointsIndices=pointsIndices,
        ν=collect(Int.(pointν)))
    return block, geometry
end

"""
    _cutCellPieces1D(params1D) -> (pieces, cellOffsets)

The 1-D integrands of one space coordinate, cell by cell, as coefficients in
s = 2y - 1 (y ∈ [0,1] the cell coordinate), multiplied by Δ (dx = Δ dy).
Same `WYYKKIntegralSymbolic` cache key as `WYYKKIntegralNumerical`.
"""
function _cutCellPieces1D(params)
    @unpack orderBspline1D, YorderBspline1Dμᶜ, YorderBspline1Dμ, μᶜs, μs, maxNode, ν,
        lᶜ_nᶜ_max, l_n_max, Δ = params
    WDerivativeOrder = Int(params.WDerivativeOrder)
    νRefPoints = params.νRefPoints
    WνOffset = params.WνOffset
    ImakeReport = params.ImakeReport
    paramsForSymbolic = @strdict orderBspline1D YorderBspline1Dμᶜ YorderBspline1Dμ μᶜs μs maxNode ν νRefPoints WνOffset WDerivativeOrder lᶜ_nᶜ_max l_n_max ImakeReport
    paramsForSymbolic["symbolic_endpoint_version"] = "segmentwise_definite_integrals_v8_overlap_partition"
    output = myProduceOrLoad(WYYKKIntegralPureSymbolic, paramsForSymbolic, "WYYKKIntegralSymbolic")
    W = output["WYYKK_integral"]
    x, Δx = W.variables[1], W.variables[2]
    nodes = [Float64(Symbolics.value(n)) for n in W.nodes]

    # cells [k, k+1] (stencil coordinates) and the segment lying in each
    segmentCell = Int[]
    for ι in 1:length(nodes)-1
        a, b = nodes[ι], nodes[ι+1]
        b - a > 1e-12 || (push!(segmentCell, typemin(Int)); continue)
        k = floor(Int, a + 1e-9)
        (isapprox(a, k; atol=1e-9) && isapprox(b, k + 1; atol=1e-9)) || throw(ArgumentError(
            "segment [$a, $b] is not a whole grid cell: sub-cell breakpoints are not supported yet"))
        push!(segmentCell, k)
    end
    cells = sort!(unique(filter(!=(typemin(Int)), segmentCell)))

    degree = max(orderBspline1D, 0) + max(YorderBspline1Dμ, 0) + max(YorderBspline1Dμᶜ, 0) +
             l_n_max + lᶜ_nᶜ_max
    sFit = [cos((2j - 1) * π / (2 * (degree + 1))) for j in 1:degree+1]   # Chebyshev nodes in s
    V = [s^p for s in sFit, p in 0:degree]
    sCheck = (-0.83, 0.11, 0.67)
    nμ, nμᶜ = length(μs), length(μᶜs)
    pieces = zeros(Float64, l_n_max + 1, lᶜ_nᶜ_max + 1, nμ, nμᶜ, length(cells), degree + 1)
    for (ι, k) in enumerate(segmentCell)
        k == typemin(Int) && continue
        iCell = searchsortedfirst(cells, k)
        xOf(s) = Δ * (k + (s + 1) / 2)
        for iμᶜ in 1:nμᶜ, iμ in 1:nμ, lc in 0:lᶜ_nᶜ_max, l in 0:l_n_max
            antiderivative = W.data[ι, 1, l+1, lc+1, iμ, iμᶜ, 1]
            integrand = Symbolics.derivative(antiderivative, x)
            f(s) = Float64(Symbolics.value(Symbolics.substitute(integrand, Dict(x => xOf(s), Δx => Δ))))
            values = [f(s) for s in sFit]
            c = V \ values
            scale = max(maximum(abs, values), 1.0)
            for s in sCheck
                abs(sum(c[p+1] * s^p for p in 0:degree) - f(s)) <= 1e-9 * scale ||
                    error("1-D integrand is not a polynomial of degree <= $degree on cell $k")
            end
            pieces[l+1, lc+1, iμ, iμᶜ, iCell, :] .= Δ .* c
        end
    end
    return pieces, cells .- ν
end

"Monomial moments (y ∈ [0,1]) → moments of s = 2y - 1, along every axis."
function _sMoments(m::AbstractArray{Float64,N}) where {N}
    p = size(m, 1) - 1
    T = [b <= a ? binomial(a, b) * 2.0^b * (-1.0)^(a - b) : 0.0 for a in 0:p, b in 0:p]
    out = copy(m)
    for d in 1:N
        out = mapslices(v -> T * v, out; dims=d)
    end
    out
end

"""
    cutCellRows(scheme, geom, nodes, models, prepared; side=:minus, model2whole=identity)

Rows of the cut nodes `nodes` (level-set grid indices) in the layout of
`prepared` (`A_unknown`, `L_known`, `R_force`), as
`Dict(row => [(column, value), ...])` for each of the three matrices.

`models` are the material arrays in the order of the equation's material
variables (e.g. `[ρ, λ, μ]`, scaled like the volume operator), defined on
the level-set grid and EXTENDED across the interface (the Taylor expansion
of a cut node still reads its whole stencil; the integral itself only covers
the material part).
"""
function cutCellRows(scheme::CutCellScheme, geom::CutCellGeometry{N}, nodes, models, prepared;
                     side::Symbol=:minus, model2whole=identity) where {N}
    N == scheme.nSpace || throw(DimensionMismatch("geometry is $N-D, scheme is $(scheme.nSpace)-D"))
    geom.maxDegree >= scheme.maxDegree || throw(ArgumentError(
        "build CutCellGeometry with maxDegree >= $(scheme.maxDegree) for this scheme"))
    own = side === :minus ? Int8(1) : Int8(0)
    p = scheme.maxDegree
    trim(m) = m[ntuple(_ -> 1:p+1, N)...]
    fullS = _sMoments(fullCellMoments(N, p))
    zeroS = zeros(Float64, ntuple(_ -> p + 1, N))
    cache = Dict{CartesianIndex{N},Array{Float64,N}}()
    function lookup(c::CartesianIndex{N})
        get!(cache, c) do
            checkbounds(Bool, geom.state, c) || return fullS          # outside the grid: as the volume
            st = geom.state[c]
            st == 2 && return _sMoments(trim(side === :minus ? geom.moments[c] :
                fullCellMoments(N, geom.maxDegree) .- geom.moments[c]))
            st == own ? fullS : zeroS
        end
    end
    _cutCellRowsWithMoments(scheme, lookup, nodes, models, prepared; model2whole=model2whole)
end

"""
    cutCellSelfCheck(scheme, prepared, models, node; model2whole=identity)

Rows of `node` computed here with every support cell uncut, compared with
the rows of the volume operator `prepared`: returns the largest relative
difference for A, L and R (should be round-off).
"""
function cutCellSelfCheck(scheme::CutCellScheme, prepared, models, node::CartesianIndex;
                          model2whole=identity)
    N = scheme.nSpace
    fullS = _sMoments(fullCellMoments(N, scheme.maxDegree))
    rows = _cutCellRowsWithMoments(scheme, _ -> fullS, [node], models, prepared; model2whole=model2whole)
    errors = map((:A, :L, :R), (prepared.A_unknown, prepared.L_known, prepared.R_force)) do key, M
        worst = 0.0
        for (row, entries) in getproperty(rows, key)
            reference = Vector(M[row, :])
            mine = zeros(length(reference))
            for (col, v) in entries
                mine[col] += v
            end
            scale = max(maximum(abs, reference), eps())
            worst = max(worst, maximum(abs, mine .- reference) / scale)
        end
        worst
    end
    (A=errors[1], L=errors[2], R=errors[3])
end

function _cutCellRowsWithMoments(scheme::CutCellScheme, lookup, nodes, models, prepared;
                                 model2whole=identity)
    N = scheme.nSpace
    T = prepared.timePointsUsedForOneStep
    scheme.pointsShape[end] == T || throw(DimensionMismatch(
        "stencil has $(scheme.pointsShape[end]) time levels, the system $T"))
    NP = prepared.NpointsSpace
    NF = prepared.NField
    NFF = prepared.NForceField
    spaceShape = prepared.spaceShape
    pointLinear = LinearIndices(spaceShape)
    modelShape = size(models[1])
    νspace = CartesianIndex(Tuple(scheme.ν[1:N]))
    stencilSpace = [CartesianIndex(Tuple(Int.(pt[1:N]))) - νspace for pt in scheme.pointsIndices]
    stencilTime = [Int(pt[end]) for pt in scheme.pointsIndices]
    rowsA = Dict{Int,Vector{Tuple{Int,Float64}}}()
    rowsL = Dict{Int,Vector{Tuple{Int,Float64}}}()
    rowsR = Dict{Int,Vector{Tuple{Int,Float64}}}()
    nodeValueCache = Dict{Any,Float64}()
    blocks = vcat([(b, false) for b in scheme.lhs], [(scheme.rhs, true)])
    entries = [_blockEntries(b, N) for (b, _) in blocks]

    for G in nodes
        lpG = pointLinear[model2whole(G)]
        materialAt = [CartesianIndex(ntuple(d -> clamp(G[d] + off[d], 1, modelShape[d]), N))
                      for off in stencilSpace]
        for (iBlock, (block, isRHS)) in enumerate(blocks)
            coefficients = _cutCellBlockCoefficients(block, scheme, lookup, G, models,
                materialAt, nodeValueCache, entries[iBlock])       # (P, NField or NForce, NExpr)
            for iExpr in 1:block.NExpr
                row = (lpG - 1) * block.NExpr + iExpr
                for p2 in eachindex(stencilSpace), iField in 1:block.NField
                    v = coefficients[p2, iField, iExpr]
                    iszero(v) && continue
                    J = G + stencilSpace[p2]
                    checkbounds(Bool, pointLinear, model2whole(J)) || continue
                    lp = pointLinear[model2whole(J)]
                    iT = stencilTime[p2]
                    if isRHS
                        col = lp + (iField - 1) * NP + (iT - 1) * NP * NFF
                        push!(get!(rowsR, row, Tuple{Int,Float64}[]), (col, v))
                    elseif iT == T
                        push!(get!(rowsA, row, Tuple{Int,Float64}[]), (lp + (iField - 1) * NP, v))
                    else
                        col = lp + (iField - 1) * NP + (iT - 1) * NP * NF
                        push!(get!(rowsL, row, Tuple{Int,Float64}[]), (col, v))
                    end
                end
            end
        end
        # rows with no entry at all in one matrix still have to be replaced
        for iExpr in 1:prepared.NField
            row = (lpG - 1) * prepared.NField + iExpr
            get!(rowsA, row, Tuple{Int,Float64}[])
            get!(rowsL, row, Tuple{Int,Float64}[])
            get!(rowsR, row, Tuple{Int,Float64}[])
        end
    end
    (A=rowsA, L=rowsL, R=rowsR)
end

"""
Node-independent bookkeeping of a block: for every α, the entries
(row of S, column of S, time factor, linear index into the spatial table)
of the contraction that `windowContraction!` performs.
"""
function _blockEntries(block::CutCellBlock, N::Int)
    nCoord = N + 1
    L = size(block.C, 2)
    Pμ, Pμᶜ = size(block.C, 3), size(block.Cc, 3)
    shapes = [size(piece)[1:4] for piece in block.spacePieces]
    slotLI = [LinearIndices(shape) for shape in shapes]
    spatialLI = LinearIndices(Tuple(prod(shape) for shape in shapes))
    entries = map(eachindex(block.alphaNodes)) do iα
        rows, cols, factors, spatial = Int[], Int[], Float64[], Int[]
        for idx in 1:size(block.table, 2)
            lᶜ = block.table[1, idx, iα]
            l = block.table[2, idx, iα]
            (lᶜ > 0 && l > 0) || continue
            kt = block.table[2+nCoord, idx, iα]
            lpt = block.table[2+2nCoord, idx, iα]
            (kt > 0 && lpt > 0) || continue
            all(d -> block.table[2+d, idx, iα] > 0 && block.table[2+nCoord+d, idx, iα] > 0, 1:N) || continue
            for iμᶜ in 1:Pμᶜ, iμ in 1:Pμ
                timeFactor = block.timeIntegral[lpt, kt, block.tableμ[nCoord, iμ], block.tableμᶜ[nCoord, iμᶜ]]
                timeFactor == 0.0 && continue
                slots = ntuple(d -> slotLI[d][block.table[2+nCoord+d, idx, iα], block.table[2+d, idx, iα],
                    block.tableμ[d, iμ], block.tableμᶜ[d, iμᶜ]], N)
                push!(rows, lᶜ + (iμᶜ - 1) * L)
                push!(cols, l + (iμ - 1) * L)
                push!(factors, timeFactor)
                push!(spatial, spatialLI[slots...])
            end
        end
        (rows=rows, cols=cols, factors=factors, spatial=spatial)
    end
    return entries
end

"Contract axis `d` of `A` (its first `k` entries) with the matrix `M` (slots × powers)."
function _contractAxis(A::AbstractArray, M::AbstractMatrix, d::Int, k::Int)
    nd = ndims(A)
    nd == 1 && return M[:, 1:k] * A[1:k]
    perm = (d, (e for e in 1:nd if e != d)...)
    B = permutedims(A, perm)
    rest = size(B)[2:end]
    Bk = reshape(collect(selectdim(B, 1, 1:k)), k, :)
    C = reshape(M[:, 1:k] * Bk, size(M, 1), rest...)
    return permutedims(C, invperm(perm))
end

"""
Spatial factor of every slot combination for one node:
T[s_1, ..., s_N] = Σ_cells Σ_a mS_cell[a] Π_d c_d(s_d, cell_d, a_d).
"""
function _spatialTable(block::CutCellBlock, moments, cellIndices)
    N = length(block.spacePieces)
    mats = [reshape(piece, prod(size(piece)[1:4]), size(piece, 5), size(piece, 6)) for piece in block.spacePieces]
    table = zeros(Float64, Tuple(size(m, 1) for m in mats))
    for K in cellIndices
        A = moments[K]
        for d in 1:N
            M = mats[d][:, K[d], :]
            A = _contractAxis(A, M, d, min(size(M, 2), size(A, d)))
        end
        table .+= A
    end
    return table
end

function _cutCellBlockCoefficients(block::CutCellBlock, scheme, lookup, G, models,
                                   materialAt, nodeValueCache, entries)
    N = scheme.nSpace
    P = size(block.C, 1)
    L = size(block.C, 2)
    Pμ, Pμᶜ = size(block.C, 3), size(block.Cc, 3)
    cellIndices = CartesianIndices(Tuple(length(o) for o in block.cellOffsets))
    moments = Dict(K => lookup(G + CartesianIndex(ntuple(d -> block.cellOffsets[d][K[d]], N)))
                   for K in cellIndices)
    spatialTable = _spatialTable(block, moments, cellIndices)
    Cmat = reshape(block.C, P, L * Pμ)
    Ccmat = reshape(block.Cc, P, L * Pμᶜ)
    coefficients = zeros(Float64, P, block.NField, block.NExpr)
    S = zeros(Float64, L * Pμᶜ, L * Pμ)
    nVars = min(length(block.vars), length(models))
    for iα in eachindex(block.alphaNodes)
        e = entries[iα]
        isempty(e.rows) && continue
        fill!(S, 0.0)
        _accumulateS!(S, e.rows, e.cols, e.factors, e.spatial, spatialTable)
        out = Ccmat * S * transpose(Cmat)            # (material point, field point)
        node = block.alphaNodes[iα]
        iField, iExpr = block.alphaField[iα], block.alphaExpr[iα]
        for p1 in 1:P
            w = _nodeValue(node, block.vars, models, materialAt[p1], nVars, nodeValueCache, iα, block)
            w == 0.0 && continue
            for p2 in 1:P
                coefficients[p2, iField, iExpr] += out[p1, p2] * w
            end
        end
    end
    return coefficients
end

# Function barrier: `spatialTable` is an N-dimensional array whose N is only known at run
# time; looping over it inside `_cutCellBlockCoefficients` dispatched dynamically on every
# entry (~9 M per node in 3-D, 2 s per node). Same arithmetic, same order.
function _accumulateS!(S::Matrix{Float64}, rows::Vector{Int}, cols::Vector{Int},
                       factors::Vector{Float64}, spatial::Vector{Int}, table::Array{Float64})
    @inbounds for k in eachindex(rows)
        S[rows[k], cols[k]] += factors[k] * table[spatial[k]]
    end
    return S
end

function _nodeValue(node, vars, models, point, nVars, cache, iα, block)
    values = ntuple(i -> Float64(models[i][point]), nVars)
    get!(cache, (objectid(block), iα, values)) do
        substitution = Dict{Any,Any}(vars[i] => values[i] for i in 1:nVars)
        Float64(Symbolics.value(Symbolics.substitute(node, substitution)))
    end
end

"""
    installRows(prepared, rows)

Replace the rows of `A_unknown`, `L_known` and `R_force` given by
`rows = (A=Dict, L=Dict, R=Dict)` (e.g. `cutCellRows`).
"""
function installRows(prepared, rows)
    A = _replace_sparse_rows(prepared.A_unknown, rows.A)
    L = _replace_sparse_rows(prepared.L_known, rows.L)
    R = _replace_sparse_rows(prepared.R_force, rows.R)
    _prepared_with_matrices(prepared, A, L, R; cut_cell_rows=sort!(collect(keys(rows.A))))
end

"""
    pinRows(prepared, nodes; model2whole=identity)

u = 0 at `nodes` (every field): identity rows, zero known/force rows. For the
nodes on the void side whose support holds no material at all.
"""
function pinRows(prepared, nodes; model2whole=identity)
    pointLinear = LinearIndices(prepared.spaceShape)
    NP, NF = prepared.NpointsSpace, prepared.NField
    rows = (A=Dict{Int,Vector{Tuple{Int,Float64}}}(), L=Dict{Int,Vector{Tuple{Int,Float64}}}(),
        R=Dict{Int,Vector{Tuple{Int,Float64}}}())
    for G in nodes
        lp = pointLinear[model2whole(G)]
        for f in 1:NF
            row = (lp - 1) * NF + f
            rows.A[row] = [(lp + (f - 1) * NP, 1.0)]
            rows.L[row] = Tuple{Int,Float64}[]
            rows.R[row] = Tuple{Int,Float64}[]
        end
    end
    installRows(prepared, rows)
end

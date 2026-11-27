# src/cutCells/aggregation.jl
#
# Small-cut-cell treatment, default chosen 2026-11-24: a node on the other
# side of the interface whose support still reaches a sliver of material is
# MERGED into a neighbour on the material side ("cell aggregation", e.g.
# Badia, Verdugo & Martín 2018):
#
#   trial side : its value is no longer an unknown but an extension of the
#                material neighbours (u_I = 2 u_J - u_K along the direction
#                I → J, or u_I = u_J when K is not available);
#   test side  : its (cut) row is added to the row of the target J, so no
#                part of the weak form is lost; its own row becomes the
#                extension constraint.
#
# The test-side merge assumes the sliver node already carries its cut-cell
# row (stage 2). With rows that are still pinned to zero (u_I = 0), call
# `applyAggregation(...; aggregate_test=false)`.

"""
    sliverAggregation(geom; side=:minus, extension=:linear, reach=1, maxSteps=reach + 1)

Void-side nodes (not on `side`) that a cut-cell row can read:

* the cut nodes on the other side (their support reaches a sliver of
  material): value extended AND row merged into the target (`aggregate=true`);
* every other node of the other side within `reach` grid steps (Chebyshev
  distance) of a cut node (the stencil half-width; 1 for OPT3): value
  extended only (`aggregate=false`) -- otherwise a cut row would read a
  pinned zero there.

The extension follows the grid direction `o` closest to the inward normal:
`J = I + k o` (k ≤ `maxSteps`, default `reach + 1`; 3-D needs more, since a
diagonal neighbour of a cut node can sit ~2.8 cells from the material) is the first node on `side`, and u_I is the Lagrange
extrapolation of degree 0 (:constant), 1 (:linear, u_I = (1+k) u_J - k u_K
with K = J + o), 2 (:quadratic) or 3 (:cubic) from J, J + o, J + 2o, ...
(lower degree where the material side has fewer nodes). Use a degree that
matches the scheme: linear for OPT3, cubic for OPT5. Returns a vector of `(node, target, extension, aggregate)`.
"""
function sliverAggregation(geom::CutCellGeometry{N}; side::Symbol=:minus,
                           extension::Symbol=:linear, reach::Int=1, maxSteps::Int=reach + 1) where {N}
    degree = Dict(:constant => 0, :linear => 1, :quadratic => 2, :cubic => 3)
    haskey(degree, extension) ||
        throw(ArgumentError("extension must be :constant, :linear, :quadratic or :cubic"))
    wanted = degree[extension]
    ls = geom.levelSet
    mask = sideMask(ls; side=side)
    cut = cutNodes(geom; side=side)
    cutSet = Set(cut)
    candidates = Set{CartesianIndex{N}}()
    for I in cut, o in CartesianIndices(ntuple(_ -> -reach:reach, N))
        J = I + o
        (checkbounds(Bool, mask, J) && !mask[J]) && push!(candidates, J)
    end
    out = NamedTuple[]
    for I in sort!(collect(candidates))
        n = levelSetNormal(ls, nodeCoordinates(ls, I))
        inward = side === :minus ? -n : n
        best = nothing
        bestK = 0
        bestScore = -Inf
        for o in CartesianIndices(ntuple(_ -> -1:1, N))
            all(iszero, Tuple(o)) && continue
            k = 0
            for s in 1:maxSteps
                J = I + s * o
                checkbounds(Bool, mask, J) || break
                if mask[J]
                    k = s
                    break
                end
            end
            k == 0 && continue
            score = dot(normalize(SVector{N,Float64}(Tuple(o) .* ls.spacing)), inward)
            if score > bestScore + 1e-12 || (abs(score - bestScore) <= 1e-12 && k < bestK)
                best, bestK, bestScore = o, k, score
            end
        end
        isnothing(best) && error("node $I has no node on side $side within $maxSteps steps (raise maxSteps)")
        J = I + bestK * best
        # Lagrange extrapolation to s = 0 (node I) from the nodes J + m·o at s = k + m,
        # m = 0..q, q = the wanted degree or less if the material side runs out
        q = 0
        while q < wanted
            K = J + (q + 1) * best
            (checkbounds(Bool, mask, K) && mask[K]) || break
            q += 1
        end
        sNodes = [Float64(bestK + m) for m in 0:q]
        weights = [(J + m * best, prod((0.0 - sNodes[j]) / (sNodes[m+1] - sNodes[j]) for j in 1:q+1 if j != m + 1; init=1.0))
                   for m in 0:q]
        push!(out, (node=I, target=J, extension=weights, aggregate=(I in cutSet)))
    end
    return out
end

"""
    applyAggregation(prepared, aggregation; model2whole=identity, aggregate_test=true)

Apply `sliverAggregation` output to a prepared linear system (the NamedTuple
of `prepareLinearSystem` / `overlapBoundaryLinearSystem`):

1. `aggregate_test`: row(target) += row(sliver node), for `A_unknown`,
   `L_known` and `R_force`, every field (only entries with `aggregate=true`);
2. row(sliver node) := u_node - Σ w u_extension (future level only;
   its `L_known` and `R_force` rows are zeroed).

`model2whole` maps the level-set grid indices to the operator grid (the
`conv.model2whole` of the operator geometry; `identity` without Cerjan
padding).
"""
function applyAggregation(prepared, aggregation; model2whole=identity,
                          aggregate_test::Bool=true)
    NField = prepared.NField
    NP = prepared.NpointsSpace
    point_linear = LinearIndices(prepared.spaceShape)
    residual = LinearIndices((NField, NP))
    nrows = size(prepared.A_unknown, 1)
    column(lp, field) = lp + (field - 1) * NP

    A, L, R = prepared.A_unknown, prepared.L_known, prepared.R_force
    if aggregate_test
        Is, Js, Vs = Int[], Int[], Float64[]
        for agg in aggregation
            (hasproperty(agg, :aggregate) && !agg.aggregate) && continue
            lp = point_linear[model2whole(agg.node)]
            lt = point_linear[model2whole(agg.target)]
            for f in 1:NField
                push!(Is, residual[f, lt]); push!(Js, residual[f, lp]); push!(Vs, 1.0)
            end
        end
        P = sparse(1:nrows, 1:nrows, ones(nrows), nrows, nrows) +
            sparse(Is, Js, Vs, nrows, nrows)
        A, L, R = P * A, P * L, P * R
    end

    constrained = Int[]
    Ic, Jc, Vc = Int[], Int[], Float64[]
    for agg in aggregation
        lp = point_linear[model2whole(agg.node)]
        for f in 1:NField
            r = residual[f, lp]
            push!(constrained, r)
            push!(Ic, r); push!(Jc, column(lp, f)); push!(Vc, 1.0)
            for (node, w) in agg.extension
                push!(Ic, r); push!(Jc, column(point_linear[model2whole(node)], f)); push!(Vc, -w)
            end
        end
    end
    keep = ones(nrows)
    keep[constrained] .= 0.0
    D = spdiagm(0 => keep)
    A = D * A + sparse(Ic, Jc, Vc, nrows, size(A, 2))
    L = D * L
    R = D * R
    dropzeros!(A); dropzeros!(L); dropzeros!(R)

    extras = (aggregation_rows=sort!(unique(constrained)), aggregation=aggregation)
    if isdefined(@__MODULE__, :_prepared_with_matrices)
        return _prepared_with_matrices(prepared, A, L, R; extras...)
    end
    return merge(prepared, (A_unknown=A, L_known=L, R_force=R), extras)
end

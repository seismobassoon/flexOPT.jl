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
    sliverAggregation(geom; side=:minus, extension=:linear)

For every node of `cutNodes(geom; side)` that is NOT on `side`, the target
neighbour `J` on `side` (the stencil neighbour whose direction is closest to
the inward normal) and the extension `[(node, weight), ...]` that defines its
value. Returns a vector of `(node, target, extension)`.
"""
function sliverAggregation(geom::CutCellGeometry{N}; side::Symbol=:minus,
                           extension::Symbol=:linear) where {N}
    extension in (:linear, :constant) ||
        throw(ArgumentError("extension must be :linear or :constant"))
    ls = geom.levelSet
    mask = sideMask(ls; side=side)
    out = NamedTuple[]
    for I in cutNodes(geom; side=side)
        mask[I] && continue
        n = levelSetNormal(ls, nodeCoordinates(ls, I))
        inward = side === :minus ? -n : n
        best = nothing
        bestScore = -Inf
        for o in CartesianIndices(ntuple(_ -> -1:1, N))
            all(iszero, Tuple(o)) && continue
            J = I + o
            (checkbounds(Bool, mask, J) && mask[J]) || continue
            direction = normalize(SVector{N,Float64}(Tuple(o) .* ls.spacing))
            score = dot(direction, inward)
            if score > bestScore
                best, bestScore = o, score
            end
        end
        isnothing(best) && error("node $I has no neighbour on side $side to merge into")
        J = I + best
        weights = [(J, 1.0)]
        if extension === :linear
            K = J + best
            if checkbounds(Bool, mask, K) && mask[K]
                weights = [(J, 2.0), (K, -1.0)]
            end
        end
        push!(out, (node=I, target=J, extension=weights))
    end
    return out
end

"""
    applyAggregation(prepared, aggregation; model2whole=identity, aggregate_test=true)

Apply `sliverAggregation` output to a prepared linear system (the NamedTuple
of `prepareLinearSystem` / `overlapBoundaryLinearSystem`):

1. `aggregate_test`: row(target) += row(sliver node), for `A_unknown`,
   `L_known` and `R_force`, every field;
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

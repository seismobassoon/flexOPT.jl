# src/cutCells/simplexIntegration.jl
#
# Exact integration of polynomials on N-simplices (N = 1, 2, 3) with the
# Grundmann & Möller (1978) rules. A rule with parameter s integrates every
# polynomial of total degree <= 2s+1 exactly (some weights are negative, which
# is harmless for the low degrees used here).
#
# Checked against ∫_triangle x^a z^b = a! b! / (a+b+2)! and the tetrahedron
# analogue (dev/cutCells/run.jl, step `geometry`).

const _GM_RULE_CACHE = Dict{Tuple{Int,Int},Tuple{Vector{Vector{Float64}},Vector{Float64}}}()

"All `parts`-tuples of non-negative integers summing to `total`."
function _compositions(total::Int, parts::Int)
    parts == 1 && return [[total]]
    out = Vector{Vector{Int}}()
    for k in 0:total
        for rest in _compositions(total - k, parts - 1)
            push!(out, vcat(k, rest))
        end
    end
    return out
end

"""
    grundmannMollerRule(n, s) -> (barycentric points, weights)

Rule of degree 2s+1 on the n-simplex. The weights are normalised so that
`n! * |T| * Σ w f(Σ λ_j v_j)` is the integral over a simplex `T` of volume `|T|`.
"""
function grundmannMollerRule(n::Int, s::Int)
    get!(_GM_RULE_CACHE, (n, s)) do
        d = 2s + 1
        points = Vector{Vector{Float64}}()
        weights = Float64[]
        for i in 0:s
            w = (-1)^i * 2.0^(-2s) * Float64(d + n - 2i)^d /
                (factorial(big(i)) * factorial(big(d + n - i)))
            for β in _compositions(s - i, n + 1)
                push!(points, [(2β[j] + 1) / (d + n - 2i) for j in 1:n+1])
                push!(weights, Float64(w))
            end
        end
        (points, weights)
    end
end

"Volume (length, area) of the simplex with vertices `V[1..N+1]`."
function simplexVolume(V::AbstractVector{SVector{N,Float64}}) where {N}
    M = [V[j+1][i] - V[1][i] for i in 1:N, j in 1:N]
    return abs(det(M)) / factorial(N)
end

"Rule parameter s exact for tensor degree `p` per axis in N dimensions."
gmParameterForTensorDegree(N::Int, p::Int) = max(0, cld(N * p - 1, 2))

"""
    addSimplexMonomialMoments!(m, V, p; sign=1.0)

`m[a_1+1, ..., a_N+1] += sign * ∫_S y_1^a_1 ... y_N^a_N dy` for every
`0 <= a_i <= p`, over the simplex `S` with vertices `V`.
"""
function addSimplexMonomialMoments!(m::AbstractArray{Float64,N},
                                    V::AbstractVector{SVector{N,Float64}},
                                    p::Int; sign::Float64=1.0) where {N}
    vol = simplexVolume(V)
    vol == 0.0 && return m
    points, weights = grundmannMollerRule(N, gmParameterForTensorDegree(N, p))
    scale = sign * factorial(N) * vol
    powers = zeros(Float64, p + 1, N)
    for (λ, w) in zip(points, weights)
        x = sum(λ[j] * V[j] for j in 1:N+1)
        for d in 1:N
            powers[1, d] = 1.0
            for k in 1:p
                powers[k+1, d] = powers[k, d] * x[d]
            end
        end
        for a in CartesianIndices(m)
            v = 1.0
            for d in 1:N
                v *= powers[a[d], d]
            end
            m[a] += scale * w * v
        end
    end
    return m
end

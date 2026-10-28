# dev/spherical1DelasticKernel/lib/sgtKernels1D.jl
#
# Frequency-domain waveform Fréchet kernels from the two SGT sides
# (Zhao & Chevrot 2011, eqs 18-19, here written in the frequency domain and
# for a transversely isotropic medium with radial symmetry axis):
#
#   δu_n(R, ω) = ∫ [ ω_c² δρ u^R·u^S  -  ε^R : δC : ε^S ] dV
#
#   u^S, ε^S : displacement and strain at Q of the source S (moment tensor)
#   u^R, ε^R : displacement and strain at Q of a unit force at R along n
#              (= by reciprocity, the receiver Green's function)
#   ω_c = ω - iω_I.
#
# (sign: -ω²ρu - ∇·(C:∇u) = f  ⇒  δu = G ⋆ (ω²δρ u + ∇·(δC:∇u)); the
# ∇·(δC:∇u) term is integrated by parts, no surface term at a free surface.)
#
# Love parameters, local basis (r, θ, φ):
#   C_rrrr = C, C_θθθθ = C_φφφφ = A, C_θθφφ = A - 2N, C_rrθθ = C_rrφφ = F,
#   C_rθrθ = C_rφrφ = L, C_θφθφ = N.
# Kernel densities K_p with δu = ∫ Σ_p K_p δp dV:
#   K_ρ = ω_c² u^R·u^S
#   K_C = -ε^R_rr ε^S_rr
#   K_A = -(ε^R_θθ + ε^R_φφ)(ε^S_θθ + ε^S_φφ)
#   K_F = -[ε^R_rr (ε^S_θθ + ε^S_φφ) + (ε^R_θθ + ε^R_φφ) ε^S_rr]
#   K_L = -4 (ε^R_rθ ε^S_rθ + ε^R_rφ ε^S_rφ)
#   K_N = -4 ε^R_θφ ε^S_θφ + 2 (ε^R_θθ ε^S_φφ + ε^R_φφ ε^S_θθ)
# Isotropic: K_λ = K_A + K_C + K_F, K_μ = 2K_A + 2K_C + K_L + K_N (λ, μ),
# K_κ = K_λ, K_μ|κ = K_μ - (2/3) K_λ.  With complex moduli (Q fixed) the
# elastic-modulus kernels are multiplied by the attenuation factor c(ω, Q).

"""
    kernelDensities(εR, εS, uR, uS, ωc) -> NamedTuple (ρ, A, C, F, L, N, λ, μ, κ, μκ)

All 3×3 strains / 3-vectors in the same local (r, θ, φ) basis at Q.
"""
function kernelDensities(εR, εS, uR, uS, ωc)
    hR = εR[2, 2] + εR[3, 3]; hS = εS[2, 2] + εS[3, 3]
    Kρ = ωc^2 * sum(uR .* uS)
    KC = -εR[1, 1] * εS[1, 1]
    KA = -hR * hS
    KF = -(εR[1, 1] * hS + hR * εS[1, 1])
    KL = -4 * (εR[1, 2] * εS[1, 2] + εR[1, 3] * εS[1, 3])
    KN = -4 * εR[2, 3] * εS[2, 3] + 2 * (εR[2, 2] * εS[3, 3] + εR[3, 3] * εS[2, 2])
    Kλ = KA + KC + KF
    Kμ = 2KA + 2KC + KL + KN
    return (ρ=Kρ, A=KA, C=KC, F=KF, L=KL, N=KN, λ=Kλ, μ=Kμ, κ=Kλ, μκ=Kμ - 2Kλ / 3)
end

"""
    waveformKernelAt(diskS, diskR, iω, S, M, Rcv, n, Q)

Kernel densities at Q (frequency index iω) for the moment tensor M at S
(components in S's local geographic basis) and the receiver component n
(unit vector in R's local basis, see `receiverDirection`).  diskS: :moment
disk at rs = S.r; diskR: :force disk at rs = Rcv.r; both must contain Q.r and
share the frequency list.
"""
function waveformKernelAt(diskS, diskR, iω, S::PointGeo, M, Rcv::PointGeo, n, Q::PointGeo; order::Integer=5)
    εS = sgtAt(diskS, iω, S, M, Q; which=:strain, order=order)
    uS = sgtAt(diskS, iω, S, M, Q; which=:disp, order=order)
    εR = sgtAt(diskR, iω, Rcv, n, Q; which=:strain, order=order)
    uR = sgtAt(diskR, iω, Rcv, n, Q; which=:disp, order=order)
    ωc = diskS.ω[iω] - im * diskS.ωI
    return kernelDensities(εR, εS, uR, uS, ωc)
end

"""
    waveformKernelAtInterp(itpS, itpR, S, M, Rcv, n, Q; side=:upper)

Same as `waveformKernelAt`, Q anywhere inside the disk radii (interpolated;
itpS/itpR = diskInterpolants of the :moment and :force disks at one ω).
"""
function waveformKernelAtInterp(itpS, itpR, S::PointGeo, M, Rcv::PointGeo, n, Q::PointGeo; side::Symbol=:upper)
    εS = sgtAtInterp(itpS, S, M, Q; which=:strain, side=side)
    uS = sgtAtInterp(itpS, S, M, Q; which=:disp, side=side)
    εR = sgtAtInterp(itpR, Rcv, n, Q; which=:strain, side=side)
    uR = sgtAtInterp(itpR, Rcv, n, Q; which=:disp, side=side)
    d = itpS.disk
    ωc = d.ω[itpS.iω] - im * d.ωI
    return kernelDensities(εR, εS, uR, uS, ωc)
end

"""
    shellBornIntegral(diskS, diskR, iω, S, M, Rcv, n; r_nodes, r_weights,
                      δ::Function, nθ, nφ)

∫ Σ_p K_p δp dV over a spherical shell with radial quadrature (r_nodes km
-- must be disk radii -- and weights), Gauss–Legendre in cos θ (nθ) and
uniform in φ (nφ).  `δ(r)` returns a NamedTuple of perturbations among
(ρ, A, C, F, L, N) (missing ones = 0).  Returns δu_n(ω).  With
`interp = (S=diskInterpolants(diskS, iω), R=diskInterpolants(diskR, iω))`
the r_nodes may be any radii inside the disks (on-the-fly interpolation).
"""
function shellBornIntegral(diskS, diskR, iω, S::PointGeo, M, Rcv::PointGeo, n;
                           r_nodes, r_weights, δ::Function, nθ::Integer, nφ::Integer, order::Integer=5,
                           interp=nothing)
    x, wx = gaussLegendre(nθ)
    φs = [2π * (j - 1) / nφ for j in 1:nφ]; wφ = 2π / nφ
    total = zero(ComplexF64)
    for (r, wr) in zip(r_nodes, r_weights)
        d = δ(r)
        acc = zero(ComplexF64)
        for (xi, wi) in zip(x, wx), φ in φs
            Q = PointGeo(r, acos(xi), φ)
            K = isnothing(interp) ? waveformKernelAt(diskS, diskR, iω, S, M, Rcv, n, Q; order=order) :
                waveformKernelAtInterp(interp.S, interp.R, S, M, Rcv, n, Q)
            v = zero(ComplexF64)
            for p in (:ρ, :A, :C, :F, :L, :N)
                hasproperty(d, p) && (v += getfield(K, p) * getfield(d, p))
            end
            acc += wi * v
        end
        total += wr * r^2 * acc * wφ
    end
    return total
end

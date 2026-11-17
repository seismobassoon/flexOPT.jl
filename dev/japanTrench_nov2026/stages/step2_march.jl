# dev/japanTrench_nov2026/stages/step2_march.jl
# Runs ONLY on a cache miss of step2_opt.jl; its last expression is the Dict saved.
#
# Predictor-corrector of src propagateLinearSystem(...; scheme=:predictor_corrector,
# predictor_mode=:lumped_future), written out as in Greenland so that only the
# rock-top nodes and two vertical sections are stored:
#   b = -L u_known + R_sheet f_sheet   (only the force columns of the sheet nodes)
#   u ← P⁻¹ b ;  repeat: u ← u + ω P⁻¹ (b - A u)

Np = preparedOPT.NpointsSpace
NF = preparedOPT.NField
Afuture = preparedOPT.A_unknown
Lknown = preparedOPT.L_known
timedPredictor = @timed prepareLumpedFuturePredictor(Afuture, preparedOPT)
lumpedPredictor = timedPredictor.value
@show extrema(lumpedPredictor.condition_estimates)

wholeLI = LinearIndices(preparedOPT.spaceShape)
padLower = CartesianIndex(Tuple(optPadding[1, :]))
toWhole(I::CartesianIndex{3}) = wholeLI[I + padLower]

# ---------------------------------------------------------------- sheet source
pw = PLANE_WAVE
pw.type in (:P, :SV, :SH) || error("PLANE_WAVE.type must be :P, :SV or :SH")
k0 = argmin(abs.(zBox .+ pw.depth))
z0 = zBox[k0]
σz = 1.5dx
layerRock = rockMask[:, :, k0]
@assert count(layerRock) > 0.9 * Nx * Ny "the sheet layer must be in rock"
vRef = 1e3 * median((pw.type === :P ? VpModel : VsModel)[:, :, k0][layerRock])   # m/s
pSlow = sind(pw.incidence) / vRef                                                   # s/m
eh = (sind(pw.azimuth), cosd(pw.azimuth))          # horizontal propagation (x = east, y = north)
px, py = pSlow * eh[1], pSlow * eh[2]
polar = pw.type === :SV ? (eh[2], -eh[1], 0.0) : (0.0, 0.0, 1.0)   # a in f = ∇m × a (unused for P)
fPeak = 1 / pw.period
t0 = 1.5 * pw.period                                # Ricker centre at the first sheet node

inner = BOX_HALF_WIDTH - pw.taper
tukey(s) = abs(s) <= inner ? 1.0 : (abs(s) >= BOX_HALF_WIDTH ? 0.0 :
    0.5 * (1 + cos(π * (abs(s) - inner) / pw.taper)))
tukeyDot(s) = (abs(s) <= inner || abs(s) >= BOX_HALF_WIDTH) ? 0.0 :
    -0.5 * π / pw.taper * sin(π * (abs(s) - inner) / pw.taper) * sign(s)

sHorizontal = [xBox[i] * eh[1] + yBox[j] * eh[2] for i in 1:Nx, j in 1:Ny]
sMin = minimum(sHorizontal)
sheetLin = Int[]; sheetW = Float64[]; sheetWx = Float64[]; sheetWy = Float64[]
sheetWz = Float64[]; sheetTau = Float64[]
for k in 1:Nz
    ζ = (zBox[k] - z0) / σz
    abs(ζ) <= 3 || continue
    G = exp(-ζ^2)
    Gz = -2ζ / σz * G
    for j in 1:Ny, i in 1:Nx
        tx, ty = tukey(xBox[i]), tukey(yBox[j])
        (tx > 0 && ty > 0) || continue
        rockMask[i, j, k] || continue
        push!(sheetLin, toWhole(CartesianIndex(i, j, k)))
        push!(sheetW, tx * ty * G)
        push!(sheetWx, tukeyDot(xBox[i]) * ty * G)
        push!(sheetWy, tx * tukeyDot(yBox[j]) * G)
        push!(sheetWz, tx * ty * Gz)
        push!(sheetTau, pSlow * (sHorizontal[i, j] - sMin))
    end
end
nSheet = length(sheetLin)
@printf("sheet: %d nodes at z = %.1f km (± %.1f km), V_ref = %.2f km/s, p = %.4f s/km, delays 0-%.1f s\n",
    nSheet, z0 / 1e3, 3σz / 1e3, vRef / 1e3, pSlow * 1e3, maximum(sheetTau))

NFP, NFF = preparedOPT.NforcePoints, preparedOPT.NForceField
sheetColumns = [sheetLin[q] + (c - 1) * NFP + (s - 1) * NFP * NFF for q in 1:nSheet, c in 1:3, s in 1:3]
Rsheet = preparedOPT.R_force[:, vec(sheetColumns)]              # node fastest, then c, then s
# force density f [N/m³] in the scaled PDE: f' = f Δt² / ρ0  [m]
sourceScale = pw.amplitude * dtOPT^2 / rho0
forceValues = zeros(nSheet, 3, 3)

# ∇m / M at one node: (w_x G r - w G p_x r', w_y G r - w G p_y r', w G_z r)
function sheetForce!(F, slot, t, W, Wx, Wy, Wz, Tau, t0, px, py, fPeak, scale, isP, polar)
    ax, ay, az = polar
    @inbounds for q in eachindex(W)
        τ = t - t0 - Tau[q]
        a = (π * fPeak * τ)^2
        e = exp(-a)
        r = (1 - 2a) * e
        rd = 2π^2 * fPeak^2 * τ * (2a - 3) * e
        gx = Wx[q] * r - W[q] * px * rd
        gy = Wy[q] * r - W[q] * py * rd
        gz = Wz[q] * r
        if isP                                          # f = -∇m
            fx, fy, fz = -gx, -gy, -gz
        else                                            # f = ∇m × a
            fx, fy, fz = gy * az - gz * ay, gz * ax - gx * az, gx * ay - gy * ax
        end
        F[q, 1, slot] = scale * fx
        F[q, 2, slot] = scale * fy
        F[q, 3, slot] = scale * fz
    end
    return F
end
sheetIsP = pw.type === :P

# ---------------------------------------------------------------- outputs
surfaceLinear = [toWhole(CartesianIndex(i, j, topRock[i, j])) for i in 1:Nx, j in 1:Ny]
ewLinear = [toWhole(CartesianIndex(i, jc, k)) for i in 1:Nx, k in 1:Nz]
nsLinear = [toWhole(CartesianIndex(ic, j, k)) for j in 1:Ny, k in 1:Nz]
nOut = optSteps ÷ traceStride + 1
nSnapMax = optSteps ÷ snapStride + 1
surfaceOut = zeros(Float32, 3, Nx, Ny, nOut)
ewOut = zeros(Float32, 3, Nx, Nz, nSnapMax)
nsOut = zeros(Float32, 3, Ny, Nz, nSnapMax)
outTimes = zeros(nOut); snapTimesOut = zeros(nSnapMax)
iOut = 1; iSnap = 1                                   # slot 1 = t = 0 (all zeros)

knownField = zeros(Np, NF, 2)            # (point, field, slot): slot 2 = newest
u = zeros(Np * NF)                       # column = point + Np (field-1)
rhs = zeros(size(Afuture, 1))
residual = similar(rhs)
correction = similar(u)
iterationsPerStep = zeros(Int, optSteps)
residualPerStep = zeros(optSteps)
stoppedEarly = false
lastStep = optSteps
progressStride = max(1, optSteps ÷ 25)
marchStart = time()

for n in 1:optSteps
    # f at t = (n-2)Δt, (n-1)Δt, nΔt for the solve of u(nΔt) (Kirishima convention)
    for s in 1:3
        sheetForce!(forceValues, s, (n - 3 + s) * dtOPT, sheetW, sheetWx, sheetWy, sheetWz,
            sheetTau, t0, px, py, fPeak, sourceScale, sheetIsP, polar)
    end
    mul!(rhs, Lknown, vec(knownField), -1.0, 0.0)
    mul!(rhs, Rsheet, vec(forceValues), 1.0, 1.0)
    applyLumpedFuturePredictor!(u, lumpedPredictor, rhs)
    rhsNorm = max(norm(rhs), eps())
    used = 0
    for it in 1:optCorrectorIterations
        mul!(residual, Afuture, u)
        residual .= rhs .- residual
        norm(residual) / rhsNorm <= optCorrectorTolerance && break
        applyLumpedFuturePredictor!(correction, lumpedPredictor, residual)
        u .+= optCorrectorRelaxation .* correction
        used = it
    end
    mul!(residual, Afuture, u)
    residual .= rhs .- residual
    iterationsPerStep[n] = used
    residualPerStep[n] = norm(residual) / rhsNorm

    amplitude = maximum(abs, u)
    if !isfinite(amplitude) || amplitude > blowupLimit
        @warn "stopping: displacement blew up" n amplitude
        global stoppedEarly = true
        global lastStep = n - 1
        break
    end
    knownField[:, :, 1] .= @view knownField[:, :, 2]
    knownField[:, :, 2] .= reshape(u, Np, NF)
    if n % traceStride == 0
        global iOut += 1
        outTimes[iOut] = n * dtOPT
        for j in 1:Ny, i in 1:Nx, c in 1:3
            surfaceOut[c, i, j, iOut] = u[surfaceLinear[i, j] + (c - 1) * Np]
        end
    end
    if n % snapStride == 0
        global iSnap += 1
        snapTimesOut[iSnap] = n * dtOPT
        for k in 1:Nz, c in 1:3
            for i in 1:Nx
                ewOut[c, i, k, iSnap] = u[ewLinear[i, k] + (c - 1) * Np]
            end
            for j in 1:Ny
                nsOut[c, j, k, iSnap] = u[nsLinear[j, k] + (c - 1) * Np]
            end
        end
    end
    if n % progressStride == 0 || n == optSteps
        elapsed = time() - marchStart
        @printf("  step %6d / %d  t = %6.1f s  max|u| = %.2e m  P-C it %d res %.1e  (%.0f s, ETA %.0f s)\n",
            n, optSteps, n * dtOPT, amplitude, used, residualPerStep[n], elapsed,
            elapsed / n * (optSteps - n))
    end
end
marchTime = time() - marchStart

Dict{String,Any}(
    "surfaceTimes" => outTimes[1:iOut], "surfaceU" => surfaceOut[:, :, :, 1:iOut],
    "snapTimes" => snapTimesOut[1:iSnap],
    "secEW" => ewOut[:, :, :, 1:iSnap], "secNS" => nsOut[:, :, :, 1:iSnap],
    "rawDt" => dtOPT,
    "sheet" => (z0=z0, sigma=σz, vRef=vRef, p=pSlow, t0=t0, nodes=nSheet,
                maxDelay=maximum(sheetTau), sMin=sMin),
    "sourceScale" => sourceScale, "rho0" => rho0,
    "pcDiagnostics" => (iterations=iterationsPerStep[1:max(lastStep, 1)],
                        residuals=residualPerStep[1:max(lastStep, 1)]),
    "stoppedEarly" => stoppedEarly,
    "timings" => merge(operatorTimings, (predictor=timedPredictor.time, march=marchTime)),
    "wholeShape" => preparedOPT.spaceShape,
)

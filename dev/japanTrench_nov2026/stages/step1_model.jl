# dev/japanTrench_nov2026/stages/step1_model.jl
# Material model on the box, cached with myProduceOrLoad (point sampling, as
# the Kirishima box):
#   getParamsAndTopo(...; velocity_model=:NIED)
#     * planet1D (PREM, by radius) everywhere: ρ, Vp, Vs, Qμ, Qκ;
#     * NIED Vp, Vs where the archive covers the point with confidence <= 0.8;
#     * GMT earth_relief: air above the topography, sea water between the
#       sea floor and sea level.
# Then (not cached, cheap): Vs of solid cells raised to VS_MIN_SOLID.

modelParams = Dict{String,Any}(
    "stage" => "model",
    "codeVersion" => CODE_VERSION["model"],
    "boxCentre" => BOX_CENTRE, "axis" => BOX_AXIS_AZIMUTH,
    "halfWidth" => BOX_HALF_WIDTH, "dx" => dx, "altMin" => altMin, "altMax" => altMax,
    "planet1D" => string(planet1D.DEFAULT_1D_MODEL[]),
    "velocity" => "NIED", "niedSource" => velocitySourceKey,
    "niedConfidenceMax" => NIED_CONFIDENCE_MAX,
    "niedOutside" => string(NIED_OUTSIDE), "niedLowConfidence" => string(NIED_LOW_CONFIDENCE),
    "topoPrecisionKm" => TOPO_PRECISION_KM, "water" => WATER, "air" => AIR,
)
devRedo!("model", modelParams, JAPANTRENCH_CACHE_DIR, "model")
modelCacheFile = devCacheFile(modelParams, JAPANTRENCH_CACHE_DIR, "model")
println("model cache: ", modelCacheFile, isfile(modelCacheFile) ? "  (hit)" : "  (miss -> computing)")
modelHash = commonBatchs.hash_parameters(modelParams)

modelResult = myProduceOrLoad(modelParams, JAPANTRENCH_CACHE_DIR, "model") do cfg
    m = getParamsAndTopo(
        geoBox, copy(box.effectiveRadii), TOPO_PRECISION_KM;   # copy: getParamsAndTopo edits the radii
        VpWater=WATER.vp, ρWater=WATER.ρ, VpAir=AIR.vp, ρAir=AIR.ρ,
        velocity_model=:NIED,
        nied_source=niedVelocitySource,
        nied_confidence_max=NIED_CONFIDENCE_MAX,
        nied_outside=NIED_OUTSIDE,
        nied_low_confidence=NIED_LOW_CONFIDENCE,
    )
    Dict{String,Any}(
        "ρ" => m.ρ, "Vp" => m.Vpv, "Vs" => m.Vsv, "Qμ" => m.Qμ, "Qκ" => m.Qκ,
        "niedMask" => Array(m.nied_mask),
        "niedPMask" => Array(m.nied_P_mask), "niedSMask" => Array(m.nied_S_mask),
        "niedPConfidence" => m.nied_P_confidence, "niedSConfidence" => m.nied_S_confidence,
    )
end
ρModel = modelResult["ρ"]
VpModel = modelResult["Vp"]
VsRaw = modelResult["Vs"]
QμModel, QκModel = modelResult["Qμ"], modelResult["Qκ"]
niedMask = BitArray(modelResult["niedMask"])
@assert size(ρModel) == (Nx, Ny, Nz)

# ---- phases
airMask = ρModel .<= AIR_DENSITY_CUTOFF
waterMask = .!airMask .& (VsRaw .<= 0)
rockMask = .!airMask .& (VsRaw .> 0)

# ---- minimum solid Vs (Yamaya et al. 2024: 1.0 km/s), Vp kept
VsModel = copy(VsRaw)
slow = rockMask .& (VsRaw .< VS_MIN_SOLID)
VsModel[slow] .= VS_MIN_SOLID
badK = rockMask .& (VpModel .^ 2 .<= (4 / 3) .* VsModel .^ 2)   # bulk modulus <= 0
@printf("Vs clamp: %d solid cells (%.2f %%) raised to %.2f km/s\n",
    count(slow), 100count(slow) / max(1, count(rockMask)), VS_MIN_SOLID)
any(badK) && @warn("$(count(badK)) solid cells have Vp² <= 4/3 Vs² after the clamp (negative bulk modulus)")

# ---- column structure: rock, then (water), then air, from the bottom up
@assert all(@view rockMask[:, :, 1]) "the bottom layer must be all rock (lower altMin)"
@assert all(@view airMask[:, :, end]) "the top layer must be all air (raise TOPO_MAX)"
topRock = zeros(Int, Nx, Ny)
columnOrderOK = trues(Nx, Ny)
for j in 1:Ny, i in 1:Nx
    k = findlast(@view rockMask[i, j, :])
    topRock[i, j] = k
    kw = findlast(@view waterMask[i, j, :])
    okRock = all(@view rockMask[i, j, 1:k])
    okWater = isnothing(kw) || all(@view waterMask[i, j, k+1:kw])
    columnOrderOK[i, j] = okRock && okWater
end
nOrderProblems = count(.!columnOrderOK)
nOrderProblems > 0 && @warn("$nOrderProblems columns are not rock / water / air from the bottom up")

# solid-top elevation (m): land surface or sea floor, at the grid's resolution
surfaceElevation = [0.5 * (altBox[i, j, topRock[i, j]] + altBox[i, j, topRock[i, j]+1]) for i in 1:Nx, j in 1:Ny]
waterThickness = dropdims(sum(waterMask; dims=3); dims=3) .* dx
niedColumnFraction = [count(@view(niedMask[i, j, :]) .& @view(rockMask[i, j, :])) /
                      max(1, count(@view rockMask[i, j, :])) for i in 1:Nx, j in 1:Ny]

vpMaximum = maximum(VpModel[rockMask])
@printf("model: %.1f %% rock, %.1f %% water, %.1f %% air; NIED on %.1f %% of the rock cells\n",
    100count(rockMask) / length(rockMask), 100count(waterMask) / length(rockMask),
    100count(airMask) / length(rockMask), 100count(niedMask .& rockMask) / count(rockMask))
@printf("       solid top from %.2f to %.2f km; water up to %.1f km (%d cells)\n",
    minimum(surfaceElevation) / 1e3, maximum(surfaceElevation) / 1e3, maximum(waterThickness) / 1e3, count(waterMask))
@printf("       rock: Vp %.2f-%.2f, Vs %.2f-%.2f km/s, ρ %.2f-%.2f g/cm³, Qμ %.0f-%.0f\n",
    extrema(VpModel[rockMask])..., extrema(VsModel[rockMask])..., extrema(ρModel[rockMask])...,
    extrema(QμModel[rockMask])...)
count(waterMask) > 0 && @warn("pure-water cells (Vs = 0) exist; μ = 0 is untested in the elastic OPT operator so far")

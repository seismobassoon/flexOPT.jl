# dev/freeSurfacePlaneWave/stages/opt_propagate.jl
# :direct marching (no predictor-corrector) -- deliberately the simplest
# possible time scheme, to keep the number of suspects small while
# checking the free-surface closure itself.
optPadding = isnothing(bcOPT.cerjan) ? zeros(Int, 2, 2) : cerjan_padding(bcOPT.cerjan)
xOPT = range(first(xOPTPhysical) - optPadding[1,1]*dxOPT; step=dxOPT, length=preparedOPT.spaceShape[1])
zOPT = range(first(zOPTPhysical) - optPadding[1,2]*dxOPT; step=dxOPT, length=preparedOPT.spaceShape[2])
sourceIZOPT = argmin(abs.(zOPTPhysical .- sourcePosition.z))
sourceLinearsOPT = [LinearIndices(preparedOPT.spaceShape)[CartesianIndex(ix + optPadding[1,1], sourceIZOPT + optPadding[1,2])]
    for ix in sourceNodesX]
optSteps = ceil(Int, duration / dtOPT)
optOutputStride = max(1, round(Int, outputSampling / dtOPT))
nPastForceLevels = preparedOPT.timePointsUsedForOneStep - 1
sourceTimesOPT = ((1 - nPastForceLevels):optSteps) .* dtOPT
waveletOPT = rickerSource.(sourceTimesOPT)
sourceScaleOPT = sourceForce * dtOPT^2 / (rho0 * dxOPT^2)   # scaled-coordinates convention
sourceOPT = zeros(Float64, preparedOPT.NforcePoints, preparedOPT.NForceField, length(sourceTimesOPT))
for linear in sourceLinearsOPT
    sourceOPT[linear, 2, :] .+= sourceScaleOPT .* waveletOPT
end
propagationOPT = propagateLinearSystem(preparedOPT, optSteps, dtOPT;
    sourceFull=sourceOPT, output_stride=optOutputStride, blowup_limit=1e8,
    solver_name="OPT3", scheme=:direct)
@assert !propagationOPT.stopped_early
uxOPT = propagationOPT.history[:, :, 1, :]
uzOPT = propagationOPT.history[:, :, 2, :]
optTimes = propagationOPT.times

sourceOPT5 = zeros(Float64, preparedOPT5.NforcePoints, preparedOPT5.NForceField, length(sourceTimesOPT))
for linear in sourceLinearsOPT
    sourceOPT5[linear, 2, :] .+= sourceScaleOPT .* waveletOPT
end
propagationOPT5 = propagateLinearSystem(preparedOPT5, optSteps, dtOPT;
    sourceFull=sourceOPT5, output_stride=optOutputStride, blowup_limit=1e8,
    solver_name="OPT5", scheme=:direct)
@assert !propagationOPT5.stopped_early
uxOPT5 = propagationOPT5.history[:, :, 1, :]
uzOPT5 = propagationOPT5.history[:, :, 2, :]
opt5Times = propagationOPT5.times
@show dtOPT size(uzOPT) maximum(abs, uzOPT) maximum(abs, uzOPT5)
# Sanity check specific to this setup: the problem is x-independent by
# construction, so every column should carry (numerically) the same trace.
midColumnOPT = cld(size(uzOPT, 1), 2)
xSpreadOPT = maximum(abs, uzOPT .- uzOPT[midColumnOPT:midColumnOPT, :, :])
xSpreadOPT5 = maximum(abs, uzOPT5 .- uzOPT5[midColumnOPT:midColumnOPT, :, :])
@show xSpreadOPT xSpreadOPT5

Dict{String,Any}(
    "xOPT" => collect(xOPT), "zOPT" => collect(zOPT), "dxOPT" => dxOPT, "dtOPT" => dtOPT,
    "uzOPT" => uzOPT, "uxOPT" => uxOPT, "optTimes" => optTimes,
    "uzOPT5" => uzOPT5, "uxOPT5" => uxOPT5, "opt5Times" => opt5Times,
)

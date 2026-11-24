# dev/freeSurfacePlaneWave/stages/opt_propagate.jl
# :direct marching (no predictor-corrector) -- deliberately the simplest
# possible time scheme, to keep the number of suspects small while
# checking the free-surface closure itself.
#
# NO FORCING AT ALL. The plane wave is injected as an EXACT initial
# condition (config.jl's planeWaveUz), evaluated at t=0 and at the past
# time levels the scheme needs, then propagated with knownForce identically
# zero throughout. propagateLinearSystem's own `initialCondition` kwarg is
# scalar-only (it just floods the whole field with one constant), so this
# is a small hand-rolled analogue of its :direct branch -- same primitives
# (flexOPT.prepareConstantMatrix, prepared.b_fun!, prepared.b_template,
# prepared.known_lhs_template / known_rhs_template, all read directly off
# `src/numSolvers/linearSolution.jl` and used here unmodified), kept here
# in dev/ only so the IC can be a full spatial field. The free-surface
# closure itself (opt_build.jl, the thing actually under test) is
# untouched -- only how energy enters the system changes.
optPadding = isnothing(bcOPT.cerjan) ? zeros(Int, 2, 2) : cerjan_padding(bcOPT.cerjan)
xOPT = range(first(xOPTPhysical) - optPadding[1,1]*dxOPT; step=dxOPT, length=preparedOPT.spaceShape[1])
zOPT = range(first(zOPTPhysical) - optPadding[1,2]*dxOPT; step=dxOPT, length=preparedOPT.spaceShape[2])
# Void above the free surface must stay EXACTLY zero displacement, matching
# the same convention FD3's own stepper enforces every step (ux_next/
# uz_next is only ever written at material points -- see
# step_elastic_wave_2d! in src/elasticWave2D.jl). The physical threshold
# z<=0 identifies the void region directly off the (padding-aware) zOPT
# coordinates, so this needs no separate material-mask plumbing: the
# bottom Cerjan padding is still solid medium (just damped), only the top
# void cap (present only in the "normal" step) is excluded.
materialRowOPT = applyFreeSurface ? (collect(zOPT) .<= 0.0) : trues(length(zOPT))
planeWaveRowOPT(iz, t) = materialRowOPT[iz] ? planeWaveUz(zOPT[iz], t) : 0.0

if !@isdefined(PLANEWAVE_EDGE_PIN_WIDTH)
    const PLANEWAVE_EDGE_PIN_WIDTH = Ref(1)
end

function propagateFreePlaneWave(prepared, Nt, dt; output_stride=1, blowup_limit=Inf, solver_name="OPT")
    edgePinWidthOPT = min(PLANEWAVE_EDGE_PIN_WIDTH[], prepared.spaceShape[1] ÷ 2)
    edgeColumnsOPT = vcat(1:edgePinWidthOPT, (prepared.spaceShape[1] - edgePinWidthOPT + 1):prepared.spaceShape[1])
    nKnownTime = prepared.timePointsUsedForOneStep - 1
    knownField = zeros(Float64, size(prepared.known_lhs_template))
    knownForce = zeros(Float64, size(prepared.known_rhs_template))   # identically zero: no forcing anywhere
    unknownField = zeros(Float64, prepared.NpointsSpace, prepared.NField)
    linear = LinearIndices(prepared.spaceShape)
    for point in CartesianIndices(prepared.spaceShape)
        iz = point[2]
        idx = linear[point]
        unknownField[idx, 2] = planeWaveRowOPT(iz, 0.0)   # field 1 = ux stays 0 (no shear, x-independent)
        for level in 1:nKnownTime
            # level=nKnownTime (most recent) -> t=0, matching unknownField
            # above; level=1 (oldest) -> t=-(nKnownTime-1)*dt. Exactly the
            # FD3 leapfrog's own (ux_previous=t-dt, ux=t) convention,
            # generalised to nKnownTime>1: the roll `knownField[:,1:end-1].=
            # knownField[:,2:end]; knownField[:,end].=unknownField` after
            # each step is what keeps this true for every later step too.
            tPast = -(nKnownTime - level) * dt
            knownField[idx, 2, level] = planeWaveRowOPT(iz, tPast)
        end
    end
    _, factor = flexOPT.prepareConstantMatrix(prepared; sparse_output=true)
    b = copy(prepared.b_template)
    snapshots = Array{Float64}[reshape(copy(unknownField), prepared.spaceShape..., prepared.NField)]
    times = Float64[0.0]
    stored_lookup = Set(unique(vcat(0, collect(output_stride:output_stride:Nt), Nt)))
    stopped_early = false
    for step in 1:Nt
        knownInputs = vcat(vec(knownField), vec(knownForce))
        prepared.b_fun!(b, knownInputs)
        unknownField .= reshape(factor \ b, prepared.NpointsSpace, prepared.NField)
        # Pin the two lateral (x) edge columns to the exact analytical
        # solution at the time just solved for -- same reasoning as
        # fd_propagate.jl's pinLateralEdgesFD!: nothing here ever imposes a
        # genuine lateral boundary condition (no periodic BC, no x-Cerjan,
        # by the original "nothing propagates laterally" design), so
        # whatever the solver's own edge treatment amounts to otherwise
        # leaks a spurious x-gradient inward once the IC is nonzero
        # everywhere. This keeps the edges exact regardless of domain
        # width, and only ever overwrites the two edge columns -- it
        # cannot mask or soften an interior/free-surface blowup.
        currentTime = step * dt
        for ix in edgeColumnsOPT, iz in 1:prepared.spaceShape[2]
            edgeIdx = linear[CartesianIndex(ix, iz)]
            unknownField[edgeIdx, 1] = 0.0
            unknownField[edgeIdx, 2] = planeWaveRowOPT(iz, currentTime)
        end
        amplitude = maximum(abs, unknownField)
        if !isfinite(amplitude) || amplitude > blowup_limit
            @warn "stopping propagation" solver_name step amplitude blowup_limit
            stopped_early = true
            break
        end
        if nKnownTime > 0
            nKnownTime > 1 && (knownField[:, :, 1:end-1] .= knownField[:, :, 2:end])
            knownField[:, :, end] .= unknownField
        end
        if step in stored_lookup
            push!(snapshots, reshape(copy(unknownField), prepared.spaceShape..., prepared.NField))
            push!(times, step * Float64(dt))
        end
    end
    history = cat(snapshots...; dims=length(prepared.spaceShape) + 2)
    (; history, times, stopped_early)
end

optSteps = ceil(Int, duration / dtOPT)
optOutputStride = max(1, round(Int, outputSampling / dtOPT))

# NOTE (2026-11-16, on request): these used to be hard
# @assert !stopped_early checks -- the "normal" (actual free-surface) step's
# OPT3 closure reliably hits blowup_limit around step ~880/1018 (the
# already-established :additive elasticTractionFree2D closure bug;
# see this benchmark's own earlier findings), which used to crash the
# whole run before FD3/SPECFEM2D/analytical ever got compared. Downgraded
# to a @warn so compare.jl still renders analytical/FD3/OPT3/OPT5/
# SPECFEM2D together for "normal" -- OPT3's truncated, divergent trace
# is itself the point being shown, not something to hide. "nosurface"
# has never hit this (no free-surface closure exercised there at all).
propagationOPT = propagateFreePlaneWave(preparedOPT, optSteps, dtOPT;
    output_stride=optOutputStride, blowup_limit=1e8, solver_name="OPT3")
propagationOPT.stopped_early && @warn "OPT3 stopped early (known free-surface-closure blowup) -- its trace is truncated/divergent in the comparison figure" step=length(propagationOPT.times)
uxOPT = propagationOPT.history[:, :, 1, :]
uzOPT = propagationOPT.history[:, :, 2, :]
optTimes = propagationOPT.times

propagationOPT5 = propagateFreePlaneWave(preparedOPT5, optSteps, dtOPT;
    output_stride=optOutputStride, blowup_limit=1e8, solver_name="OPT5")
propagationOPT5.stopped_early && @warn "OPT5 stopped early (known free-surface-closure blowup) -- its trace is truncated/divergent in the comparison figure" step=length(propagationOPT5.times)
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

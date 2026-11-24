# dev/freeSurfacePointSource/seismicWavePCEnv/runSeismicWavePC.jl -- SeismicWavePC.jl
# run of the point-force/flat-free-surface benchmark, in an ISOLATED Julia
# environment (its own Project.toml/Manifest.toml under this directory, NOT
# flexOPT's shared one at the repo root).
#
# Why isolated: SeismicWavePC.jl pulls in GeophysicalModelGenerator, which
# pins Interpolations/DIVAnd/ChainRulesCore/Tullio versions that conflict
# with flexOPT's own pinned Makie/GLMakie and its own use of Tullio in
# src/flexOPT.jl (confirmed: job 037, "Unsatisfiable requirements... resolved
# by compatibility requirements with Makie... GLMakie..."). Keeping
# SeismicWavePC in its own environment sidesteps that entirely -- this
# script is never `include`d from, and never Pkg.activate()s, flexOPT's own
# Project.toml.
#
# This reproduces the physical setup of
# notebooks/HomogeneousElastic2DBenchmark_FreeSurface.ipynb cell 3 (quickRun
# values only) / dev/freeSurfacePointSource/run.jl: vp0=6000, vs0=3464,
# rho0=2700 (m/s, m/s, kg/m3), sourcePosition=(x=-10km,z=-20km),
# dx=500m, domainHalfWidth=60km, duration=9s, receiverZ=-10km line at
# receiverX=[-30,-20,-10,0,10] km, Ricker peak frequency
# 1.5*vs0/(16*500m) ~ 0.6495 Hz. A faithful rebuild of dev/freeSurfacePointSource/smoke.jl's
# already-worked-out SeismicWavePC API usage (ForceSource/Station/simulate),
# just now in its own environment so Pkg.add actually succeeds.
#
# Three conventions do NOT translate one-to-one from the other solvers (see
# application/forBoris/README.md for the first two):
#   1. Force amplitude: OpenSWPC's own internal M0 normalization, not a
#      direct Newton import of sourceForce=1e10 N/m. -> traces are saved
#      RAW (not pre-normalized); the comparison script normalizes every
#      solver's trace to its own peak per station already, so this doesn't
#      bias the shape/timing comparison.
#   2. Source time function: SeismicWavePC's built-in :ricker (trise
#      convention) vs this benchmark's own rickerSource formula -- close but
#      not byte-identical.
#   3. Source geometry: a true point force here (SeismicWavePC's
#      ForceSource), vs the other solvers' Gaussian-smoothed source
#      (sigma=500m, cutoff=2sigma). Left as-is for this first pass.
#
# Run (note --project points at THIS directory, not the repo root):
#   julia --project=dev/freeSurfacePointSource/seismicWavePCEnv -t auto \
#       dev/freeSurfacePointSource/seismicWavePCEnv/runSeismicWavePC.jl

import Pkg
Pkg.activate(@__DIR__)

neededRegistryPackages = ["GeophysicalModelGenerator", "JLD2", "KernelAbstractions"]
for pkg in neededRegistryPackages
    isnothing(Base.find_package(pkg)) && Pkg.add(pkg)
end
if isnothing(Base.find_package("SeismicWavePC"))
    println("Adding SeismicWavePC.jl (github.com/JuliaGeodynamics/SeismicWavePC.jl)...")
    try
        Pkg.add("SeismicWavePC")
    catch
        Pkg.add(Pkg.PackageSpec(url="https://github.com/JuliaGeodynamics/SeismicWavePC.jl.git"))
    end
end
Pkg.instantiate()

using SeismicWavePC, GeophysicalModelGenerator, JLD2
using KernelAbstractions: CPU

# ---- physical setup (quickRun values, identical to run.jl / smoke.jl) ----
vp0, vs0, rho0 = 6000.0, 3464.0, 2700.0   # m/s, m/s, kg/m3
sourcePosition = (x=-10e3, z=-20e3)       # m
dx = 500.0                                 # m
domainHalfWidth = 60e3                     # m
duration = 9.0                             # s
receiverZ = -10e3
receiverX = [-30e3, -20e3, -10e3, 0.0, 10e3]
referencePointsPerSWavelength = 16
rickerSharpness = 1.5
sourceFrequency = rickerSharpness * vs0 / (referencePointsPerSWavelength * dx)
@show sourceFrequency

# ---- model: homogeneous half-space, air above z=0, grid+fields in km, km/s ----
xkm = collect((sourcePosition.x - domainHalfWidth):dx:(sourcePosition.x + domainHalfWidth)) ./ 1e3
zkm = collect((sourcePosition.z - domainHalfWidth):dx:(sourcePosition.z + domainHalfWidth)) ./ 1e3
model = CartData(xyz_grid(xkm, 0, zkm))
X, Y, Z = xyz_grid(xkm, 0, zkm)
isAir = Z .> 0.0
rhoField = ifelse.(isAir, 1.2, rho0)
vpField = ifelse.(isAir, 0.0, vp0 / 1e3)   # km/s
vsField = ifelse.(isAir, 0.0, vs0 / 1e3)   # km/s
model = addfield(model, (; rho=rhoField, vp=vpField, vs=vsField))

trise = 1.0 / sourceFrequency
source = ForceSource(sourcePosition.x / 1e3, 0.0, sourcePosition.z / 1e3; fz=1.0, t0=0.0, trise=trise)
stations = [Station("S$(i)", xr / 1e3, 0.0, receiverZ / 1e3) for (i, xr) in enumerate(receiverX)]

dt = 0.3 * (dx / 1e3) / (vp0 / 1e3 * sqrt(2))
nt = ceil(Int, duration / dt)
@show dt nt size(model)

@time result = simulate(model; dt, nt, sources=source, stf=:ricker, stations, backend=CPU(), verbose=true)
display(result)
@show fieldnames(typeof(result))
@show result.components
@show size(result.u)
@show result.stations

# ---- extract uz traces at the receiver line, defensively (field layout TBD
#      ahead of first actual run -- print everything needed to debug if the
#      shape guesses below are wrong) ----
compIndex = findfirst(c -> c in (:z, :uz, :Z, "z", "uz", "Z"), result.components)
isnothing(compIndex) && error("Could not find a vertical component in result.components = $(result.components)")

nStations, nComponents, nSamples = length(stations), length(result.components), length(result.t)
uDims = size(result.u)
@show uDims nStations nComponents nSamples

timeAxis = Float64.(result.t)
traces = Vector{NamedTuple{(:time, :values),Tuple{Vector{Float64},Vector{Float64}}}}(undef, length(receiverX))
for i in eachindex(receiverX)
    values = if uDims == (nStations, nComponents, nSamples)
        Float64.(result.u[i, compIndex, :])
    elseif uDims == (nSamples, nComponents, nStations)
        Float64.(result.u[:, compIndex, i])
    elseif uDims == (nComponents, nStations, nSamples)
        Float64.(result.u[compIndex, i, :])
    elseif uDims == (nComponents, nSamples, nStations)
        Float64.(result.u[compIndex, :, i])
    elseif uDims == (nStations, nSamples, nComponents)
        Float64.(result.u[i, :, compIndex])
    elseif uDims == (nSamples, nStations, nComponents)
        Float64.(result.u[:, i, compIndex])
    else
        error("Unrecognized result.u shape $(uDims) for (nStations=$nStations, nComponents=$nComponents, nSamples=$nSamples); inspect result.u manually")
    end
    traces[i] = (time=timeAxis, values=values)
end

outdir = joinpath(@__DIR__, "..", "outputs", "seismicWavePC")
mkpath(outdir)
jldsave(joinpath(outdir, "seismicwavepc_traces.jld2");
    vp0, vs0, rho0, sourcePosition, sourceFrequency, dx, domainHalfWidth, duration,
    receiverX, receiverZ, traces, dt, nt)
println("saved -> ", joinpath(outdir, "seismicwavepc_traces.jld2"))
println("DONE")

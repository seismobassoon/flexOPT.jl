# dev/cutCells/seismicWavePCCases.jl -- SeismicWavePC.jl (OpenSWPC-type FD, staggered grid)
# runs of the topography cases, in the ISOLATED environment of
# dev/freeSurfacePointSource/seismicWavePCEnv (see runSeismicWavePC.jl there for why).
#
# SeismicWavePC has no cut cells: its free surface is the air/solid contact of a
# staircase model (air: vp = vs = 0, rho = 1.2 where z > z_s(x)), i.e. it is a ZIGZAG
# code. Same grid spacing (500 m), same medium, same receivers as flexOPT/SPECFEM2D;
# force along n̂ = (-sin θ, cos θ).
#
# Known convention differences (unchanged from the flat benchmark):
#   - force amplitude: OpenSWPC's internal normalisation, not 1e10 N/m;
#   - STF: SeismicWavePC's built-in :ricker (trise = 1/f0), not gcRicker;
#   - true point force, not the 500 m Gaussian distribution.
# They are fixed ONCE, globally, by compareAll.jl: one factor and one lag measured on
# tilt0 against SPECFEM2D and applied unchanged to every case (stated on the figures).
# No per-case or per-trace normalisation.
#
# Run:
#   julia --project=dev/freeSurfacePointSource/seismicWavePCEnv -t auto \
#       dev/cutCells/seismicWavePCCases.jl [tilt0 tilt10 tilt20 tilt30 hill]
# Output: dev/cutCells/outputs/references/seismicwavepc_<case>.jld2

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "freeSurfacePointSource", "seismicWavePCEnv"))
include(joinpath(@__DIR__, "geometryCases.jl"))
using SeismicWavePC, GeophysicalModelGenerator, JLD2, Printf
using KernelAbstractions: CPU

const SWPC_DX = 500.0

function componentTrace(result, comp, i)
    nS, nC, nT = length(result.stations), length(result.components), length(result.t)
    d = size(result.u)
    d == (nS, nC, nT) && return Float64.(result.u[i, comp, :])
    d == (nT, nC, nS) && return Float64.(result.u[:, comp, i])
    d == (nC, nS, nT) && return Float64.(result.u[comp, i, :])
    d == (nC, nT, nS) && return Float64.(result.u[comp, :, i])
    d == (nS, nT, nC) && return Float64.(result.u[i, :, comp])
    d == (nT, nS, nC) && return Float64.(result.u[:, i, comp])
    error("unrecognised result.u shape $d")
end

function runSeismicWavePCCase(case)
    hasFreeSurface(case) || error("SeismicWavePC runs only for cases with a free surface")
    dx = SWPC_DX
    x, z = caseGrid(case, dx)
    z = collect(first(z):dx:(last(z) + 4dx))      # a few more air layers on top
    s = caseSurface(case)
    xkm, zkm = x ./ 1e3, z ./ 1e3
    X, Y, Z = xyz_grid(xkm, 0, zkm)
    isAir = (Z .* 1e3) .> s.(X .* 1e3)
    model = CartData(xyz_grid(xkm, 0, zkm))
    model = addfield(model, (; rho=ifelse.(isAir, 1.2, GC.rho),
        vp=ifelse.(isAir, 0.0, GC.vp / 1e3), vs=ifelse.(isAir, 0.0, GC.vs / 1e3)))
    sx, sz = caseSource(case)
    n = caseForce(case)
    source = ForceSource(sx / 1e3, 0.0, sz / 1e3; fx=n[1], fz=n[2], t0=0.0, trise=1 / gcFrequency())
    receivers = caseReceivers(case)
    points = vcat(receivers[:deep], receivers[:shallow])
    stations = [Station("S$(i)", p[1] / 1e3, 0.0, p[2] / 1e3) for (i, p) in enumerate(points)]
    dt = 0.3 * (dx / 1e3) / (GC.vp / 1e3 * sqrt(2))
    nt = ceil(Int, GC.duration / dt)
    @printf("  SeismicWavePC %s: %d × %d nodes, %d air, dt = %.4f s, nt = %d\n",
        case, length(x), length(z), count(isAir), dt, nt)
    t0 = time()
    result = simulate(model; dt, nt, sources=source, stf=:ricker, stations, backend=CPU(), verbose=false)
    @printf("  SeismicWavePC %s done in %.0f s; components %s, u %s\n",
        case, time() - t0, string(result.components), string(size(result.u)))
    cx = findfirst(c -> c in (:x, :ux, :X, "x", "ux", "X"), result.components)
    cz = findfirst(c -> c in (:z, :uz, :Z, "z", "uz", "Z"), result.components)
    (isnothing(cx) || isnothing(cz)) && error("components x/z not found in $(result.components)")
    traces = Dict{String,Any}()
    nPerLine = length(GC.receiversX)
    for (j, line) in enumerate(("deep", "shallow")), ir in 1:nPerLine
        k = (j - 1) * nPerLine + ir
        ut, un = caseRotate(case, componentTrace(result, cx, k), componentTrace(result, cz, k))
        traces["$(line)_$(ir)"] = (x_km=GC.receiversX[ir] / 1e3, line=line, ut=Float64.(ut), un=Float64.(un))
    end
    outDir = joinpath(@__DIR__, "outputs", "references")
    mkpath(outDir)
    jldsave(joinpath(outDir, "seismicwavepc_$(case).jld2"); case, times=Float64.(result.t), traces,
        note="raw SeismicWavePC output (point force, built-in ricker); global calibration in compareAll.jl")
    println("  traces -> ", joinpath(outDir, "seismicwavepc_$(case).jld2"))
end

swpcCaseList = isempty(ARGS) ? ["tilt0", "tilt10", "tilt20", "tilt30", "hill"] : ARGS
for case in swpcCaseList
    println("\n== SeismicWavePC: $case")
    runSeismicWavePCCase(case)
end

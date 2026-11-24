# dev/freeSurfacePointSource/smoke.jl -- FIRST-PASS smoke test, not yet the
# cached run.jl/stages/ pipeline used elsewhere in dev/ (see dev/README.md).
# Goal: confirm SeismicWavePC.jl (github.com/JuliaGeodynamics/SeismicWavePC.jl)
# runs on the *exact* physical setup of the existing point-force/flat-free-
# surface benchmark (notebooks/HomogeneousElastic2DBenchmark_FreeSurface.ipynb,
# cell 3: vp0=6000, vs0=3464, rho0=2700 m/s & kg/m3, sourcePosition=(-10km,
# -20km), dx=500m, duration=9s, receiverZ=-10km line at
# receiverX=[-30,-20,-10,0,10] km), and produces a first overlay against the
# SPECFEM2D reference traces already cached on disk at the same 5 stations
# (data/specfem2d_benchmarks/homogeneous_flat_free_surface_dx500m_nrec23,
# stations S0016/S0022/S0017/S0023/S0018 -- confirmed via DATA/STATIONS).
#
# Two things are DELIBERATELY NOT solved here yet (flagged, not hidden):
#   1. Force amplitude: SeismicWavePC's ForceSource magnitude follows
#      OpenSWPC's own internal source normalization (see simulate()'s "M0"
#      bookkeeping in src/simulation.jl), not a direct import of our
#      sourceForce=1e10 N/m line-force convention. Traces below are each
#      normalized to their own peak for a shape/timing comparison only.
#   2. Source time function: using SeismicWavePC's built-in :ricker (peak
#      frequency 1/trise) rather than reproducing our own exact rickerSource
#      formula via a SampledSTF -- close in character, not byte-identical.
# If this smoke test looks physically sane (right arrival order/timing,
# plausible P/S/Rayleigh structure), a second pass builds the proper
# dev/freeSurfacePointSource/run.jl + stages/ pipeline with OPT3/OPT5/FD3
# re-run alongside, and tackles amplitude calibration.

import Pkg
flexopt_root = dirname(dirname(@__DIR__))
isfile(joinpath(flexopt_root, "src", "flexOPT.jl")) ||
    error("expected to run as dev/freeSurfacePointSource/smoke.jl from the repo root")
Pkg.activate(flexopt_root)

if isnothing(Base.find_package("SeismicWavePC"))
    println("Adding SeismicWavePC.jl (github.com/JuliaGeodynamics/SeismicWavePC.jl) to the project...")
    try
        Pkg.add("SeismicWavePC")
    catch
        Pkg.add(Pkg.PackageSpec(url="https://github.com/JuliaGeodynamics/SeismicWavePC.jl.git"))
    end
end

using SeismicWavePC, GeophysicalModelGenerator
using CairoMakie
using DrWatson, JLD2

include(joinpath(flexopt_root, "src", "specfemBenchmark.jl"))   # read-only: find_specfem2d_traces / read_specfem2d_trace
using .specfemBenchmark

# ---- physical setup, copied from the notebook's cell 3 (quickRun values) ----
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

println("sourceFrequency = ", sourceFrequency, " Hz")

# ---- SeismicWavePC model: homogeneous half-space, air above z=0, in km ----
xkm = collect((sourcePosition.x - domainHalfWidth):dx:(sourcePosition.x + domainHalfWidth)) ./ 1e3
zkm = collect((sourcePosition.z - domainHalfWidth):dx:(sourcePosition.z + domainHalfWidth)) ./ 1e3
model = CartData(xyz_grid(xkm, 0, zkm))
X, Y, Z = xyz_grid(xkm, 0, zkm)
isAir = Z .> 0.0
rhoField = ifelse.(isAir, 1.2, rho0)
vpField = ifelse.(isAir, 0.0, vp0 / 1e3)
vsField = ifelse.(isAir, 0.0, vs0 / 1e3)
model = addfield(model, (; rho=rhoField, vp=vpField, vs=vsField))

trise = 1.0 / sourceFrequency
source = ForceSource(sourcePosition.x / 1e3, 0.0, sourcePosition.z / 1e3; fz=1.0, t0=0.0, trise=trise)

stations = [Station("S$(i)", xr / 1e3, 0.0, receiverZ / 1e3) for (i, xr) in enumerate(receiverX)]

dt = 0.3 * (dx / 1e3) / (vp0 / 1e3 * sqrt(2))
nt = ceil(Int, duration / dt)
println("dt = ", dt, " s, nt = ", nt, " (grid ", size(model), ")")

@time result = simulate(model; dt, nt, sources=source, stf=:ricker, stations, backend=CPU(), verbose=true)
display(result)

outdir = joinpath(@__DIR__, "outputs")
mkpath(outdir)
jldsave(joinpath(outdir, "smoke_seismicwavepc.jld2"); result, sourceFrequency, dt, nt, receiverX, receiverZ)

# ---- SPECFEM2D reference traces at the same 5 stations (already cached) ----
specfemOutput = joinpath(flexopt_root, "data", "specfem2d_benchmarks",
    "homogeneous_flat_free_surface_dx500m_nrec23", "OUTPUT_FILES")
specfemStationNames = ["S0016", "S0022", "S0017", "S0023", "S0018"]  # x = -30,-20,-10,0,10 km @ z=-10km
specfemFiles = find_specfem2d_traces(specfemOutput; component=:z, network="FX", stations=specfemStationNames)
@assert length(specfemFiles) == 5
specfemTraces = [read_specfem2d_trace(f) for f in specfemFiles]

fig = Figure(size=(1000, 900))
for (i, xr) in enumerate(receiverX)
    ax = Axis(fig[i, 1]; ylabel="x=$(round(Int, xr/1e3)) km",
        xlabel=(i == length(receiverX) ? "time (s)" : ""))
    uz = result.u[:, 2, i]
    lines!(ax, result.t, uz ./ maximum(abs, uz); label="SeismicWavePC (normalized)")
    st = specfemTraces[i]
    lines!(ax, st.time, st.values ./ maximum(abs, st.values); label="SPECFEM2D velocity (normalized)", linestyle=:dash)
    i == 1 && axislegend(ax; position=:rt)
end
Label(fig[0, 1], "Point-force flat-free-surface smoke test: SeismicWavePC vs SPECFEM2D (shape/timing only, NOT calibrated amplitude)"; fontsize=14)
save(joinpath(outdir, "smoke_seismicwavepc_vs_specfem.png"), fig)
println("figure -> ", joinpath(outdir, "smoke_seismicwavepc_vs_specfem.png"))
println("DONE")

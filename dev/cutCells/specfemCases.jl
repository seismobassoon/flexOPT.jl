# dev/cutCells/specfemCases.jl -- SPECFEM2D reference runs of the topography cases
#
# Same problem as flexOPT (geometryCases.jl): homogeneous medium, the same Gaussian
# distribution of point forces on the 500 m grid nodes, the same Ricker source time
# function (external STF), force along n̂, receivers at the same points; the mesh
# follows the free surface exactly (top interface = z_s(x)); the bottom interface is
# parallel to the surface, 80 km (measured normal to a tilted surface) below it, so
# every column has the same element height; Stacey absorbing bottom/left/right.
# Displacement seismograms (seismotype 1), both components.
#
# Run (flexOPT environment, needs the SPECFEM2D installation used by
# src/specfemBenchmark.jl):
#   julia --project=. dev/cutCells/specfemCases.jl                 # tilt0 tilt10 tilt20 tilt30 hill
#   julia --project=. dev/cutCells/specfemCases.jl tilt20 hill
#   FLEXOPT_REDO_SPECFEM=1 julia --project=. dev/cutCells/specfemCases.jl ...
# Cases live in data/specfem2d_benchmarks/cutcells_<case>_dx500m; traces are collected
# in dev/cutCells/outputs/references/specfem_<case>.jld2.

if !@isdefined(flexopt_root)
    import Pkg
    global flexopt_root = abspath(joinpath(@__DIR__, "..", ".."))
    Pkg.activate(flexopt_root)
end
include(joinpath(@__DIR__, "geometryCases.jl"))
!@isdefined(specfemBenchmark) && include(joinpath(flexopt_root, "src", "specfemBenchmark.jl"))
using .specfemBenchmark, JLD2, Printf

const SPECFEM_DX = 500.0
const SPECFEM_ELEMENT = 4SPECFEM_DX          # 2 km spectral elements (5 GLL points), as the flat benchmark
specfemDt() = 0.20 * SPECFEM_DX / (sqrt(2) * GC.vp * 2)   # = flexOPT Δt (0.00589 s), as the flat benchmark

function writeParallelInterfaces(path, x, surface, depth, nzElements)
    open(path, "w") do io
        println(io, "2")
        println(io, length(x))
        for (xi, si) in zip(x, surface)
            println(io, xi, " ", si - depth)
        end
        println(io, length(x))
        for (xi, si) in zip(x, surface)
            println(io, xi, " ", si)
        end
        println(io, nzElements)
    end
    path
end

function runSpecfemCase(case; redo=get(ENV, "FLEXOPT_REDO_SPECFEM", "0") == "1")
    hasFreeSurface(case) || error("SPECFEM2D reference only for cases with a free surface")
    dx = SPECFEM_DX
    xGrid, zGrid = caseGrid(case, dx)
    s = caseSurface(case)
    θ = deg2rad(caseTheta(case))
    depth = 80e3 / cos(θ)                    # vertical thickness; 80 km normal to a tilted plane
    x = xGrid
    surface = s.(x)
    zmin = floor((minimum(surface) - depth - dx) / dx) * dx
    zmax = ceil((maximum(surface) + dx) / dx) * dx
    z = collect(zmin:dx:zmax)
    vp = fill(GC.vp, length(x), length(z))
    vs = fill(GC.vs, length(x), length(z))
    rho = fill(GC.rho, length(x), length(z))
    nxElements = cld(length(x) - 1, 4)
    nzElements = round(Int, depth / SPECFEM_ELEMENT)
    weights = caseSourceWeights(case, xGrid, zGrid)
    sources = [(x=xGrid[i], z=zGrid[k], weight=w) for (i, k, w) in weights]
    receivers = caseReceivers(case)
    points = vcat([(x=p[1], z=p[2]) for p in receivers[:deep]], [(x=p[1], z=p[2]) for p in receivers[:shallow]])
    caseDir = joinpath(flexopt_root, "data", "specfem2d_benchmarks", "cutcells_$(case)_dx$(round(Int, dx))m")
    dt = specfemDt()
    if redo || !isfile(joinpath(caseDir, "solver.log"))
        prepared = prepare_specfem2d_case(caseDir, x, z, vp, vs, rho, surface;
            sources=sources, receiver_points=points, record_at_surface_same_vertical=false,
            duration=GC.duration, dt=dt, f0=gcFrequency(), source_factor=GC.force,
            source_angle=caseTheta(case),       # SPECFEM force: (-sin a, cos a) = n̂ for a = θ
            source_type=1, seismotype=1, source_time_function=gcRicker,
            nx_elements=nxElements, nz_elements=nzElements, free_surface=true)
        writeParallelInterfaces(prepared.interfaces, x, surface, depth, nzElements)
        println("  SPECFEM2D $case: $(nxElements) × $(nzElements) elements, running ...")
        run = run_specfem2d_case(caseDir)
        @printf("  SPECFEM2D %s done in %.0f s\n", case, run.total_wall_time_s)
    else
        println("  SPECFEM2D $case: existing run reused ($caseDir)")
    end
    out = joinpath(caseDir, "OUTPUT_FILES")
    names = ["S" * lpad(i, 4, '0') for i in 1:length(points)]
    read(component) = begin
        files = find_specfem2d_traces(out; component=component, seismotype=1, network="FX", stations=names)
        byName = Dict(split(basename(f), '.')[2] => f for f in files)
        [read_specfem2d_trace(byName[n]; time_shift=1.2 / gcFrequency()) for n in names]
    end
    tx, tz = read(:x), read(:z)
    traces = Dict{String,Any}()
    nPerLine = length(GC.receiversX)
    for (j, line) in enumerate(("deep", "shallow")), ir in 1:nPerLine
        k = (j - 1) * nPerLine + ir
        ut, un = caseRotate(case, tx[k].values, tz[k].values)
        traces["$(line)_$(ir)"] = (x_km=GC.receiversX[ir] / 1e3, line=line, ut=Float64.(ut), un=Float64.(un))
    end
    times = Float64.(tz[1].time)
    outDir = joinpath(@__DIR__, "outputs", "references")
    mkpath(outDir)
    jldsave(joinpath(outDir, "specfem_$(case).jld2"); case, times, traces)
    # sanity check of the force direction: above the source u_t vanishes by symmetry (tilted plane)
    if caseKind(case) == :tilt
        tr = traces["deep_3"]
        @printf("  check: max|u_t| / max|u_n| above the source = %.3f (≈ 0 if the force angle convention is right)\n",
            maximum(abs, tr.ut) / maximum(abs, tr.un))
    end
    println("  traces -> ", joinpath(outDir, "specfem_$(case).jld2"))
end

specfemCaseList = isempty(ARGS) ? ["tilt0", "tilt10", "tilt20", "tilt30", "hill"] : ARGS
for case in specfemCaseList
    println("\n== SPECFEM2D reference: $case")
    runSpecfemCase(case)
end

# dev/cutCells/specfem3DCases.jl -- SPECFEM3D_Cartesian reference runs of the 3-D cases
#
# Same problem as flexOPT (geometryCases3D.jl): homogeneous medium, the same Gaussian
# distribution of point forces on the 500 m nodes (one FORCESOLUTION entry per node,
# factor F·w), the same Ricker (external STF, SPECFEM3D's time origin is 0 with an external
# file), force along n, receivers at the same points (USE_SOURCES_RECEIVERS_Z: z in metres).
# Internal mesher (xmeshfem3D): one layer between a flat bottom and the top interface,
# which IS the free surface (the plane, sampled on a 2-km grid, bilinear = exact for a plane).
# Stacey absorbing conditions on the 4 sides and the bottom. Serial build (NPROC = 1).
#
# SPECFEM3D: ENV SPECFEM3D_ROOT, default ../othersPackages/specfem3d next to flexOPT; needs
# bin/xmeshfem3D, bin/xgenerate_databases, bin/xspecfem3D (./configure FC=gfortran CC=gcc
# --without-mpi && make xmeshfem3D xgenerate_databases xspecfem3D).
#
# Run:  julia --project=. dev/cutCells/specfem3DCases.jl [flat3D oblique10 oblique20]
#       FLEXOPT_REDO_SPECFEM=1 ... to redo existing runs
# Cases in data/specfem3d_benchmarks/cutcells3d_<case>; traces in
# dev/cutCells/outputs/references/specfem3d_<case>.jld2 (u1, u2, un per receiver).

if !@isdefined(flexopt_root)
    import Pkg
    global flexopt_root = abspath(joinpath(@__DIR__, "..", ".."))
    Pkg.activate(flexopt_root)
end
include(joinpath(@__DIR__, "geometryCases3D.jl"))
using JLD2, Printf

const SF3 = (element=2000.0, lateral=25e3, below=25e3, dt=0.0025)
specfem3DRoot() = get(ENV, "SPECFEM3D_ROOT", normpath(joinpath(flexopt_root, "..", "othersPackages", "specfem3d")))

"Set `KEY = value` lines of a SPECFEM parameter file (keys must exist)."
function setParameters!(text::String, values::Dict{String,String})
    for (key, value) in values
        pattern = Regex("(?m)^" * key * "\\s*=.*\$")
        occursin(pattern, text) || error("parameter $key not found in template")
        text = replace(text, pattern => rpad(key, 32) * "= " * value)
    end
    text
end
fbool(b) = b ? ".true." : ".false."
fd(x) = @sprintf("%.6fd0", x)

function prepareSpecfem3DCase(case, dir)
    root = specfem3DRoot()
    isfile(joinpath(root, "bin", "xspecfem3D")) || error("SPECFEM3D not built in $root (set SPECFEM3D_ROOT)")
    s = case3DSource(case)
    snap(v) = round(v / GC3.dx) * GC3.dx
    cx, cy, cz = snap(s[1]), snap(s[2]), snap(s[3])
    L = SF3.lateral
    xmin, xmax, ymin, ymax = cx - L, cx + L, cy - L, cy + L
    zbottom = cz - SF3.below
    zs = case3DSurface(case)
    nex = round(Int, 2L / SF3.element)
    thickness = maximum(zs(x, y) for x in (xmin, xmax), y in (ymin, ymax)) - zbottom
    nz = ceil(Int, thickness / SF3.element)
    nstep = ceil(Int, GC3.duration / SF3.dt)

    data = joinpath(dir, "DATA")
    mesh = joinpath(data, "meshfem3D_files")
    mkpath(mesh)
    mkpath(joinpath(dir, "OUTPUT_FILES"))
    mkpath(joinpath(dir, "DATABASES_MPI"))
    # top interface = free surface, on a 2-km grid covering the box
    nip = nex + 1
    open(joinpath(mesh, "interface_top.dat"), "w") do io
        for j in 0:nex, i in 0:nex                  # x fastest (meshfem3D reads ix inside iy)
            @printf(io, "%.6f\n", zs(xmin + i * SF3.element, ymin + j * SF3.element))
        end
    end
    open(joinpath(mesh, "interfaces.dat"), "w") do io
        println(io, "# number of interfaces\n 1")
        println(io, "# SUPPRESS_UTM_PROJECTION NXI NETA LONG_MIN LAT_MIN SPACING_XI SPACING_ETA")
        @printf(io, " .true. %d %d %s %s %s %s\n", nip, nip, fd(xmin), fd(ymin), fd(SF3.element), fd(SF3.element))
        println(io, " interface_top.dat")
        println(io, "# number of spectral elements in the vertical direction, layer 1")
        println(io, " ", nz)
    end
    meshPar = """
    LATITUDE_MIN                    = $(fd(ymin))
    LATITUDE_MAX                    = $(fd(ymax))
    LONGITUDE_MIN                   = $(fd(xmin))
    LONGITUDE_MAX                   = $(fd(xmax))
    DEPTH_BLOCK_KM                  = $(fd(-zbottom / 1e3))
    UTM_PROJECTION_ZONE             = 11
    SUPPRESS_UTM_PROJECTION         = .true.
    INTERFACES_FILE                 = interfaces.dat
    CAVITY_FILE                     = no_cavity.dat
    NEX_XI                          = $nex
    NEX_ETA                         = $nex
    NPROC_XI                        = 1
    NPROC_ETA                       = 1
    USE_REGULAR_MESH                = .true.
    NDOUBLINGS                      = 0
    NZ_DOUBLING_1                   = 0
    NZ_DOUBLING_2                   = 0
    CREATE_ABAQUS_FILES             = .false.
    CREATE_DX_FILES                 = .false.
    CREATE_VTK_FILES                = .false.
    SAVE_MESH_AS_CUBIT              = .false.
    LOCAL_PATH                      = ./DATABASES_MPI
    THICKNESS_OF_X_PML              = 0.d0
    THICKNESS_OF_Y_PML              = 0.d0
    THICKNESS_OF_Z_PML              = 0.d0
    NMATERIALS                      = 1
    1  $(GC3.rho)  $(GC3.vp)  $(GC3.vs)  9999.  9999.  0  2
    NREGIONS                        = 1
    1  $nex  1  $nex  1  $nz  1
    """
    write(joinpath(mesh, "Mesh_Par_file"), meshPar)

    par = read(joinpath(root, "DATA", "Par_file"), String)
    par = setParameters!(par, Dict(
        "SUPPRESS_UTM_PROJECTION" => ".true.", "NPROC" => "1", "NSTEP" => string(nstep),
        "DT" => @sprintf("%.6f", SF3.dt), "MODEL" => "default", "TOPOGRAPHY" => ".false.",
        "ATTENUATION" => ".false.", "PML_CONDITIONS" => ".false.", "STACEY_ABSORBING_CONDITIONS" => ".true.",
        "STACEY_INSTEAD_OF_FREE_SURFACE" => ".false.", "USE_SOURCES_RECEIVERS_Z" => ".true.",
        "USE_FORCE_POINT_SOURCE" => ".true.", "USE_RICKER_TIME_FUNCTION" => ".false.",
        "USE_EXTERNAL_SOURCE_FILE" => ".true.", "SAVE_SEISMOGRAMS_DISPLACEMENT" => ".true.",
        "SAVE_SEISMOGRAMS_VELOCITY" => ".false.", "SAVE_SEISMOGRAMS_ACCELERATION" => ".false.",
        "MOVIE_SURFACE" => ".false.", "MOVIE_VOLUME" => ".false.", "NTSTEP_BETWEEN_OUTPUT_INFO" => "500",
        "LOCAL_PATH" => "./DATABASES_MPI"))
    write(joinpath(data, "Par_file"), par)

    # source time function sampled at (it - 1) DT, the time SPECFEM3D applies it at
    open(joinpath(data, "stf_ricker.txt"), "w") do io
        for it in 1:nstep
            @printf(io, "%.10e\n", gc3Ricker((it - 1) * SF3.dt))
        end
    end
    # one force per node of the 500 m grid inside the Gaussian support (same as flexOPT)
    g = case3DGrid(case)
    weights = case3DSourceWeights(case, g.x, g.y, g.z)
    n = case3DForce(case)
    open(joinpath(data, "FORCESOLUTION"), "w") do io
        for (isrc, (i, j, k, w)) in enumerate(weights)
            println(io, "FORCE  ", lpad(isrc, 3, '0'))
            println(io, "time shift:     0.0000")
            println(io, "f0:             $(gc3Frequency())")
            println(io, "latorUTM:       $(g.y[j])")
            println(io, "longorUTM:      $(g.x[i])")
            println(io, "depth:          $(g.z[k])")
            println(io, "source time function:            0")
            println(io, "factor force source:             $(GC3.force * w)")
            println(io, "component dir vect source E:     $(n[1])")
            println(io, "component dir vect source N:     $(n[2])")
            println(io, "component dir vect source Z_UP:  $(n[3])")
            println(io, joinpath(data, "stf_ricker.txt"))
        end
    end
    receivers = case3DReceivers(case)
    open(joinpath(data, "STATIONS"), "w") do io
        for line in (:deep, :shallow), (ir, p) in enumerate(receivers[line])
            # STA NET lat(y) lon(x) elevation burial(= z with USE_SOURCES_RECEIVERS_Z)
            @printf(io, "%s%02d FX %.4f %.4f 0.0 %.4f\n", line == :deep ? "D" : "S", ir, p.r[2], p.r[1], p.r[3])
        end
    end
    (nex=nex, nz=nz, nstep=nstep, nsources=length(weights), root=root)
end

function runSpecfem3DCase(case; redo=get(ENV, "FLEXOPT_REDO_SPECFEM", "0") == "1")
    dir = joinpath(flexopt_root, "data", "specfem3d_benchmarks", "cutcells3d_$(case)")
    if redo || !isfile(joinpath(dir, "OUTPUT_FILES", "output_solver.txt")) ||
       !occursin("End of the simulation", read(joinpath(dir, "OUTPUT_FILES", "output_solver.txt"), String))
        info = prepareSpecfem3DCase(case, dir)
        @printf("  SPECFEM3D %s: %d × %d × %d elements, %d forces, %d steps\n", case, info.nex, info.nex, info.nz, info.nsources, info.nstep)
        t = @elapsed for exe in ("xmeshfem3D", "xgenerate_databases", "xspecfem3D")
            println("    running ", exe)
            cd(dir) do
                run(pipeline(`$(joinpath(info.root, "bin", exe))`; stdout=joinpath(dir, "OUTPUT_FILES", "$(exe).log"),
                    stderr=joinpath(dir, "OUTPUT_FILES", "$(exe).err")))
            end
        end
        @printf("  SPECFEM3D %s done in %.0f s\n", case, t)
    else
        println("  SPECFEM3D $case: existing run reused ($dir)")
    end
    out = joinpath(dir, "OUTPUT_FILES")
    readTrace(name, comp) = begin
        hits = filter(f -> occursin(Regex("^FX\\.$(name)\\..X$(comp)\\.semd\$"), f), readdir(out))   # band code depends on DT
        length(hits) == 1 || error("seismogram FX.$(name).?X$(comp).semd not found (or ambiguous) in $out")
        file = joinpath(out, only(hits))
        a = [parse.(Float64, split(l)) for l in eachline(file) if !isempty(strip(l))]
        (first.(a), last.(a))
    end
    receivers = case3DReceivers(case)
    traces = Dict{String,Any}()
    times = Float64[]
    for line in (:deep, :shallow), (ir, p) in enumerate(receivers[line])
        name = @sprintf("%s%02d", line == :deep ? "D" : "S", ir)
        tx, ux = readTrace(name, "X"); _, uy = readTrace(name, "Y"); _, uz = readTrace(name, "Z")
        times = tx
        u1, u2, un = case3DRotate(case, ux, uy, uz)
        traces["$(line)_$(ir)"] = (xp_km=p.xp / 1e3, yp_km=p.yp / 1e3, line=String(line), u1=u1, u2=u2, un=un)
    end
    outDir = joinpath(@__DIR__, "outputs", "references")
    mkpath(outDir)
    jldsave(joinpath(outDir, "specfem3d_$(case).jld2"); case, times, traces)
    tr = traces["deep_3"]           # above the source: u_1 = u_2 = 0 by symmetry
    @printf("  check above the source: max(|u_1|,|u_2|)/max|u_n| = %.3f (≈ 0)\n",
        max(maximum(abs, tr.u1), maximum(abs, tr.u2)) / maximum(abs, tr.un))
    println("  traces -> ", joinpath(outDir, "specfem3d_$(case).jld2"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    for case in (isempty(ARGS) ? collect(CASES3D) : ARGS)
        println("\n== SPECFEM3D reference: $case")
        runSpecfem3DCase(case)
    end
end

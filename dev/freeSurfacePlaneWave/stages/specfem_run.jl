# dev/freeSurfacePlaneWave/stages/specfem_run.jl
# Runs ONLY on a cache miss of specfem.jl; last expression is the Dict
# saved. Builds a case directory, runs SPECFEM2D's own mesher+solver on it
# (via run_specfem2d_case, which shells out to xmeshfem2D/xspecfem2D on
# whatever machine this Julia process is running on), reads the receiver
# traces, and integrates velocity -> displacement for comparison with the
# FD3/OPT3/OPT5/analytical uz traces in compare.jl.

caseDirectory = joinpath(flexopt_root, "data", "specfem2d_benchmarks",
    "planewave_$(DEV_STEP[])")

# SPECFEM needs its own proper free-surface mesh: "normal" uses only the
# z<=0 (solid) part, with the actual flat free surface at z=0 as the mesh's
# top interface; "nosurface" keeps the FULL margin on both sides (fully
# solid, no boundary at z=0 at all), matching FD3/OPT3's own "nosurface"
# grid exactly (see config.jl).
solidZ = applyFreeSurface ? findall(z .<= 0.0) : collect(eachindex(z))
zSPECFEM = z[solidZ]
vpField = fill(vp0, nx, length(zSPECFEM))
vsField = fill(vs0, nx, length(zSPECFEM))
rhoField = fill(rho0, nx, length(zSPECFEM))
# The actual top of THIS mesh (0.0 for "normal", z0+belowMargin for
# "nosurface") -- never a hardcoded 0.0, so this is correct either way.
meshTop = last(zSPECFEM)

dtSPECFEM = dtOPT
specfemNxElements = max(4, round(Int, (last(x) - first(x)) / (4referenceDx)))
specfemNzElements = max(4, round(Int, (meshTop - first(zSPECFEM)) / (4referenceDx)))

# source=(x=0,z=z0): under initialfield, SPECFEM2D reinterprets xs/zs as
# the reference point the analytic plane-wave snapshot is built around at
# t=0 -- NOT a point/sheet source location (there is no point or sheet
# source in this test at all). zs=meshTop=0.0 (the first guess, based on a
# manual note about zs sitting on the mesh's top interface) put nearly all
# of the t=0 field within a few km of the surface and left the z=-28km
# (below z0, direct-wave-only) receiver essentially at zero -- i.e. that
# note is about the REFLECTION-timing calculation specifically, not the
# spatial reference for where the pulse itself sits at t=0. z0 is OUR OWN
# reference depth (config.jl's planeWaveUz is centered there too), so use
# it for both steps -- the one choice that is actually symmetric with our
# own closed form rather than an SPECFEM2D-manual guess.
specfemCase = prepare_specfem2d_case(
    caseDirectory, x, zSPECFEM, vpField, vsField, rhoField,
    fill(meshTop, nx);
    source=(x=0.0, z=z0),
    receiver_points=[(x=0.0, z=depth) for depth in observerDepths],
    duration=duration, dt=dtSPECFEM, f0=sourceFrequency,
    source_angle=0.0,                   # normal incidence
    # NOT passing source_time_function: that forces time_function_type=8
    # (an external STF array), which produced the right t=0 displacement
    # norm (~1.0, matching planewaveAmplitude) but still blew up by step 5
    # (max norm ~1.8e12, SIGILL) -- a type-8 array apparently isn't enough
    # for SPECFEM2D to build a Newmark-consistent initial VELOCITY field
    # under initialfield, only displacement. Falling back to the
    # template's own default time_function_type (1 = Ricker, same
    # closed-form shape as our own rickerSource, just not byte-for-byte
    # the same t0 convention) lets SPECFEM2D build displacement+velocity
    # consistently itself; f0 below still ties its frequency to ours.
    # prepare_specfem2d_case's own source_factor default (1.0e10) is sized
    # for a conventional point force/moment tensor (whose physical units
    # need huge amplification to produce an O(1) displacement). Under
    # initialfield the SOURCE file's `factor` still multiplies the
    # injected plane-wave field itself, and our own planewaveAmplitude
    # convention is already O(1) in displacement units.
    source_factor=planewaveAmplitude,
    nx_elements=specfemNxElements, nz_elements=specfemNzElements,
    free_surface=applyFreeSurface,
    output_wavefield_dumps=false,
)

# prepare_specfem2d_case only knows the conventional point/moment-tensor
# source_type values (1=force, 2=moment tensor; it validates against
# exactly those two, so source_type=4 can't be passed in directly). Patch
# the two lines SPECFEM2D itself reinterprets once initialfield is on, by
# hand, after prepare_specfem2d_case has written a normally-valid case:
#   Par_file: initialfield = .false. -> .true.
#   SOURCE:   source_type  = 1 (its own "force" meaning) -> 1 (initialfield
#             meaning: P wave + free-surface reflection/conversion) for
#             "normal", or 4 (initialfield meaning: incident P wave only --
#             no reflected/converted wave, matching the FD3/OPT3
#             "nosurface" control exactly) for "nosurface".
function _patch_parameter_line(text, key, newvalue)
    pattern = Regex("(?m)^(\\s*" * key * "\\s*=\\s*)\\S+(.*)\$")
    m = match(pattern, text)
    isnothing(m) && error("could not find a $(key) line to patch")
    text[1:prevind(text, m.offset)] * m.captures[1] * string(newvalue) *
        m.captures[2] * text[(m.offset + ncodeunits(m.match)):end]
end

parText = read(specfemCase.par_file, String)
parText = _patch_parameter_line(parText, "initialfield", ".true.")
# PML + initialfield is a known-bad combination (PML's split-field
# auxiliary memory variables start at zero and are never consistent
# with a globally-nonzero plane-wave field right at the PML boundary,
# unlike a point source whose energy genuinely travels INTO the PML
# from a quiet interior): the template's own PML_BOUNDARY_CONDITIONS=
# .true. default produced the exact same ~1.8e12 blowup by time step 5
# regardless of the source-time-function fix above. Switch to plain
# Stacey absorbing conditions instead, SPECFEM2D's own standard pairing
# for initialfield (see e.g. its plane-wave example Par_files).
parText = _patch_parameter_line(parText, "PML_BOUNDARY_CONDITIONS", ".false.")
parText = _patch_parameter_line(parText, "STACEY_ABSORBING_CONDITIONS", ".true.")
write(specfemCase.par_file, parText)

sourceType = applyFreeSurface ? 1 : 4
sourceText = read(specfemCase.source_file, String)
sourceText = _patch_parameter_line(sourceText, "source_type", sourceType)
write(specfemCase.source_file, sourceText)

if PLANEWAVE_RUN_SPECFEM2D[]
    specfemRun = run_specfem2d_case(specfemCase.case_directory)
else
    specfemRun = (output=specfemCase.output,)
end

# SPECFEM2D's own seismogram files start at t=-1.2/f0 (confirmed directly
# off the raw .semv files: first sample at -1.8475750577...s for
# f0=0.6495Hz, exactly 1.2/0.6495) -- a standard Ricker-source convention
# (room for the pulse's causal ramp-up before its own nominal t=0),
# independent of initialfield specifically. specfemCase.time_axis_shift
# is only ever set when source_time_function is passed (not the case
# here, since v3 switched to the template's own built-in Ricker), so it
# is always 0.0 on this path and this shift was never being applied --
# which is why every trace looked like it arrived ~1.85s too early.
verticalFiles = find_specfem2d_traces(specfemRun.output; component=:z)
@assert length(verticalFiles) == length(observerDepths) "expected $(length(observerDepths)) traces, found $(length(verticalFiles)) in $(specfemRun.output)"
# NOT 1.2/f0: that only re-anchors the FILE's own start-of-recording
# label (SPECFEM2D's standard "give the Ricker's causal ramp-up room"
# convention) and got the arrival ~sourceDelay-1.2/f0 (~0.46s) too early
# everywhere -- because SPECFEM2D's own analytic plane wave uses a bare,
# zero-centered Ricker (peak exactly when t=(z-zs)/vp), while our own
# rickerSource is centered at sourceDelay, not 0. Matching
# RickerBare(t_file-(z-zs)/vp) to rickerSource(t_physical-(z-z0)/vp0) =
# RickerBare([t_physical-(z-z0)/vp0]-sourceDelay) with zs=z0 gives, for
# every z, t_physical = t_file + sourceDelay -- a single global shift,
# not 1.2/f0 at all (that number only ever described the file's own
# internal start-of-recording bookkeeping, not our pulse convention).
specfemTimeShift = sourceDelay
specfemRawTracesZ = [read_specfem2d_trace(file;
    time_shift=specfemTimeShift) for file in verticalFiles]

function integrate_velocity_trace(trace)
    displacement = zeros(Float64, length(trace.values))
    displacement[2:end] .= cumsum(
        ((trace.values[1:end-1] .+ trace.values[2:end]) ./ 2) .*
        diff(trace.time))
    (time=Float64.(trace.time), values=displacement)
end
specfemTracesZ = integrate_velocity_trace.(specfemRawTracesZ)

solverLogText = read(joinpath(specfemCase.case_directory, "solver.log"), String)
specfemCFLMatch = match(r"Max CFL stability condition[^=]*=\s*([0-9.Ee+-]+)",
    solverLogText)
specfemCFL = isnothing(specfemCFLMatch) ? missing :
    parse(Float64, specfemCFLMatch.captures[1])
@show specfemCase.case_directory specfemCFL specfemCase.time_axis_shift meshTop sourceType

Dict{String,Any}(
    "specfemTracesZ" => specfemTracesZ,
    "specfemCase" => (case_directory=specfemCase.case_directory,
        time_axis_shift=specfemCase.time_axis_shift,
        nx_elements=specfemCase.nx_elements, nz_elements=specfemCase.nz_elements,
        output=specfemRun.output),
    "specfemCFL" => specfemCFL,
    "specfemObserverDepths" => collect(observerDepths),
)

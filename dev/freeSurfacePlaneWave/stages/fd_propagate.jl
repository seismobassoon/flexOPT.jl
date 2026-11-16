# dev/freeSurfacePlaneWave/stages/fd_propagate.jl
fdConfig = ElasticThreePointConfig2D(pointsInSpace=3, pointsInTime=3, supplementaryOrder=2, cfl=0.38)
model = (ρ=fill(rho0, nx, nz) ./ 1e3, Vpv=fill(vp0, nx, nz) ./ 1e3, Vsv=fill(vs0, nx, nz) ./ 1e3)
fd = prepare_elastic_wave_2d(model, (dx, dx); material_mask=material, boundary_conditions=bcFD, config=fdConfig)
fdCoordinates = elastic_wave_coordinates(x, z, fd)
fdPadding = fd.padding
# NO FORCING AT ALL. Seed the two leapfrog time levels (t=-fd.dt and t=0)
# directly from the exact closed-form plane-wave solution (config.jl's
# planeWaveUz), masked to fd.material only: void nodes stay EXACTLY zero,
# the same convention step_elastic_wave_2d! itself enforces every step
# (ux_next/uz_next is only ever written at material points -- a nonzero
# void seed would just be an artifact of the t=0/-dt snapshot, inconsistent
# with every later step and with FD3's already-SPECFEM-matching behaviour).
fdZ2D = repeat(reshape(collect(fdCoordinates.z), 1, :), size(fd.ux, 1), 1)
fd.ux_previous .= 0; fd.uz_previous .= 0
fd.ux .= 0; fd.uz .= 0
fd.uz_previous[fd.material] .= Float32.(planeWaveUz.(fdZ2D[fd.material], -fd.dt))
fd.uz[fd.material] .= Float32.(planeWaveUz.(fdZ2D[fd.material], 0.0))
fdSteps = ceil(Int, duration / fd.dt)
fdOutputStride = max(1, round(Int, outputSampling / fd.dt))
# Pin the two lateral (x) edge columns to the exact analytical solution
# after every step. step_elastic_wave_2d! only ever updates INTERIOR
# columns (its loops run i in 2:nx-1), so the edge columns are otherwise
# simply frozen forever at whatever they were last set to. That is
# invisible for a field that starts and stays at rest (the old point-
# source design, edges always exactly 0), but becomes a real (and wrong)
# boundary condition the instant the IC is nonzero everywhere, as it now
# is: the frozen edge disagrees more and more with its evolving interior
# neighbour, leaking a spurious x-gradient one column inward per step
# (this is what xSpreadFD/xSpreadOPT below are catching). Re-pinning to
# the exact solution after every step keeps the edges honest regardless
# of domain width -- nxHalfWidth does not need to be grown to outrun this
# at the wave speed.
function pinLateralEdgesFD!()
    for i in (1, size(fd.ux, 1)), k in 1:size(fd.ux, 2)
        fd.ux[i, k] = 0
        fd.uz[i, k] = fd.material[i, k] ? Float32(planeWaveUz(fdZ2D[i, k], fd.time)) : 0f0
    end
end
fdFramesX = Matrix{Float32}[copy(fd.ux)]
fdFramesZ = Matrix{Float32}[copy(fd.uz)]
fdTimes = Float64[0.0]
for step in 1:fdSteps
    step_elastic_wave_2d!(fd)
    pinLateralEdgesFD!()
    if step % fdOutputStride == 0 || step == fdSteps
        push!(fdFramesX, copy(fd.ux)); push!(fdFramesZ, copy(fd.uz)); push!(fdTimes, fd.time)
    end
end
uxFD = cat(fdFramesX...; dims=3)
uzFD = cat(fdFramesZ...; dims=3)
@show fd.dt size(uzFD) maximum(abs, uzFD)
midColumnFD = cld(size(uzFD, 1), 2)
xSpreadFD = maximum(abs, uzFD .- uzFD[midColumnFD:midColumnFD, :, :])
@show xSpreadFD
Dict{String,Any}("uxFD"=>uxFD, "uzFD"=>uzFD, "fdTimes"=>fdTimes,
    "fdX"=>collect(fdCoordinates.x), "fdZ"=>collect(fdCoordinates.z),
    "fdPadding"=>fdPadding, "fdDt"=>Float64(fd.dt))

# dev/freeSurfacePlaneWave/stages/fd_propagate.jl
fdConfig = ElasticThreePointConfig2D(pointsInSpace=3, pointsInTime=3, supplementaryOrder=2, cfl=0.38)
model = (ρ=fill(rho0, nx, nz) ./ 1e3, Vpv=fill(vp0, nx, nz) ./ 1e3, Vsv=fill(vs0, nx, nz) ./ 1e3)
fd = prepare_elastic_wave_2d(model, (dx, dx); material_mask=material, boundary_conditions=bcFD, config=fdConfig)
fdCoordinates = elastic_wave_coordinates(x, z, fd)
fdPadding = fd.padding
sourcePhysicalZ = argmin(abs.(z .- sourcePosition.z))
sourceRowFD = sourcePhysicalZ + fdPadding[1,2]
sourceNodesFD = [CartesianIndex(ix + fdPadding[1,1], sourceRowFD) for ix in 1:nx]
fdSteps = ceil(Int, duration / fd.dt)
fdOutputStride = max(1, round(Int, outputSampling / fd.dt))
fdFramesX = Matrix{Float32}[copy(fd.ux)]
fdFramesZ = Matrix{Float32}[copy(fd.uz)]
fdTimes = Float64[0.0]
for step in 1:fdSteps
    for index in sourceNodesFD
        elasticWave2D.add_ricker_source!(fd, index; f0=sourceFrequency, t0=sourceDelay,
            amplitude=sourceForce, component=:z, source_kind=:force)
    end
    step_elastic_wave_2d!(fd)
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

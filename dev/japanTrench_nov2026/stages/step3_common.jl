# dev/japanTrench_nov2026/stages/step3_common.jl
# Coordinates, fence geometry and colour scales shared by step3 (figures) and step4 (movie).

xKm, yKm, zKm = xBox ./ 1e3, yBox ./ 1e3, zBox ./ 1e3
zTopKm = [zKm[topRock[i, j]] for i in 1:Nx, j in 1:Ny]
XS = [xKm[i] for i in 1:Nx, j in 1:Ny]; YS = [yKm[j] for i in 1:Nx, j in 1:Ny]
XE = [xKm[i] for i in 1:Nx, k in 1:Nz]; YE = fill(yKm[jc], Nx, Nz); ZE = [zKm[k] for i in 1:Nx, k in 1:Nz]
XN = fill(xKm[ic], Ny, Nz); YN = [yKm[j] for j in 1:Ny, k in 1:Nz]; ZN = [zKm[k] for j in 1:Ny, k in 1:Nz]
airEW = .!materialFS[:, jc, :]
airNS = .!materialFS[ic, :, :]
componentName = ("u_x (≈ east)", "u_y (≈ north)", "u_z (up)")
surfaceIndexAt(t) = argmin(abs.(surfaceTimes .- t))

function fenceColours(c, isnap)
    E = Float64.(secEW[c, :, :, isnap]); E[airEW] .= NaN
    N = Float64.(secNS[c, :, :, isnap]); N[airNS] .= NaN
    S = Float64.(surfaceU[c, :, :, surfaceIndexAt(snapTimes[isnap])])
    return S, E, N
end

"symmetric colour limit: 99th percentile of |u_c| over sections and surface"
function robustLimit(c)
    v = sort!(abs.(vcat(vec(secEW[c, :, :, :]), vec(secNS[c, :, :, :]), vec(surfaceU[c, :, :, :]))))
    return max(v[clamp(round(Int, 0.99 * length(v)), 1, length(v))], eps(Float32))
end

fenceAxis(pos; title="") = Axis3(pos; aspect=(1, 1, 0.45), azimuth=-0.62π, elevation=0.20π,
    xlabel="x, east (km)", ylabel="y, north (km)", zlabel="z (km)", title=title)

"E-W and N-S sections + rock-top surface; colours may be Observables"
function drawFence!(ax, S, E, N; clim)
    kw = (colormap=:balance, colorrange=(-clim, clim), nan_color=:transparent, shading=NoShading)
    surface!(ax, XE, YE, ZE; color=E, kw...)
    surface!(ax, XN, YN, ZN; color=N, kw...)
    h = surface!(ax, XS, YS, zTopKm; color=S, kw...)
    lines!(ax, xKm, fill(yKm[jc], Nx), fill(sheet.z0 / 1e3, Nx); color=:black, linestyle=:dash, linewidth=1)
    lines!(ax, fill(xKm[ic], Ny), yKm, fill(sheet.z0 / 1e3, Ny); color=:black, linestyle=:dash, linewidth=1)
    return h
end

waveLabel = @sprintf("%s plane wave, incidence %g°, towards azimuth %g°, Ricker %g s; ocean load %s; %s (dx = %.1f km)",
    PLANE_WAVE.type, PLANE_WAVE.incidence, PLANE_WAVE.azimuth, PLANE_WAVE.period, OCEAN_LOAD,
    japanTrenchCase, dx / 1e3)

"wiggle record section: traces (ntraces, nt) at `offsets` (km), common `scale` (km per m)"
function wiggles!(ax, times, offsets, traces, scale)
    for (o, tr) in zip(offsets, eachrow(traces))
        lines!(ax, times, o .+ scale .* tr; color=:black, linewidth=0.6)
    end
    return ax
end

# dev/japanTrench_nov2026/stages/step3_figures.jl
# outputs/step3_<case>/:
#   ocean_load          bathymetry (GMT) and the loaded top-cell density
#   fence_uz, fence_ux  3-D fences (E-W + N-S sections + rock-top surface) at 6 times
#   record_EW, record_NS wiggle record sections (x, y, z) on the rock top along y = 0 / x = 0,
#                       with the dashed moveout of the imposed horizontal slowness p

# ---------------------------------------------------------------- ocean load
figLoad = Figure(size=(1150, 520))
axB = Axis(figLoad[1, 1]; aspect=DataAspect(), xlabel="x, east (km)", ylabel="y, north (km)",
    title="GMT bathymetry / topography")
hB = heatmap!(axB, xKm, yKm, bathymetry ./ 1e3; colormap=:oleron, colorrange=(-8, 8))
contour!(axB, xKm, yKm, bathymetry; levels=[0.0], color=:black, linewidth=1.2)
Colorbar(figLoad[2, 1], hB; vertical=false, label="km")
axL = Axis(figLoad[1, 2]; aspect=DataAspect(), xlabel="x, east (km)", ylabel="y, north (km)",
    title="ρ of the top rock cell, with the ocean load")
hL = heatmap!(axL, xKm, yKm, [ρLoaded[i, j, topRock[i, j]] for i in 1:Nx, j in 1:Ny]; colormap=:viridis)
contour!(axL, xKm, yKm, bathymetry; levels=[0.0], color=:white, linewidth=1.2)
Colorbar(figLoad[2, 2], hL; vertical=false, label="g/cm³")
devSave(figLoad, "ocean_load")

# ---------------------------------------------------------------- fences
firstUseful = findfirst(>=(sheet.t0), snapTimes)
firstUseful = isnothing(firstUseful) ? 1 : firstUseful
fenceFrames = unique(round.(Int, range(firstUseful, length(snapTimes); length=6)))
for c in (3, 1)
    clim = robustLimit(c)
    fig = Figure(size=(1500, 950))
    local h
    for (n, isnap) in enumerate(fenceFrames)
        ax = fenceAxis(fig[fld1(n, 3), mod1(n, 3)]; title=@sprintf("t = %.0f s", snapTimes[isnap]))
        h = drawFence!(ax, fenceColours(c, isnap)...; clim=clim)
    end
    Colorbar(fig[1:2, 4], h; label=componentName[c] * " (m), clipped at the 99th percentile")
    Label(fig[0, :], waveLabel; fontsize=14)
    devSave(fig, c == 3 ? "fence_uz" : "fence_ux")
end

# ---------------------------------------------------------------- record sections
eh = (sind(PLANE_WAVE.azimuth), cosd(PLANE_WAVE.azimuth))
function recordSection(along)
    idx = 1:2:(along == :x ? Nx : Ny)                       # every other node
    offsets = along == :x ? xKm[idx] : yKm[idx]
    tr(c) = along == :x ? Float64.(surfaceU[c, idx, jc, :]) : Float64.(surfaceU[c, ic, idx, :])
    amax = maximum(maximum(abs, tr(c)) for c in 1:3)
    scale = 1.5 * 2dx / 1e3 / max(amax, eps())
    # moveout of the imposed slowness through the peak |u_z| of the centre trace
    centre = along == :x ? ic : jc
    uzCentre = along == :x ? surfaceU[3, ic, jc, :] : surfaceU[3, ic, jc, :]
    tc = surfaceTimes[argmax(abs.(uzCentre))]
    sOf(o) = along == :x ? o * eh[1] : o * eh[2]
    fig = Figure(size=(1500, 800))
    for c in 1:3
        ax = Axis(fig[1, c]; xlabel="time (s)", ylabel=(along == :x ? "x, east (km)" : "y, north (km)"),
            title=componentName[c])
        wiggles!(ax, surfaceTimes, offsets, tr(c), scale)
        lines!(ax, tc .+ sheet.p * 1e3 .* (sOf.(offsets) .- sOf(0.0)), offsets;
            color=:red, linestyle=:dash, linewidth=1)
        xlims!(ax, 0, last(surfaceTimes))
    end
    Label(fig[0, :], (along == :x ? "E-W line through the centre (y = 0)" : "N-S line through the centre (x = 0)") *
        @sprintf(" — common scale (max %.2e m); red: imposed p = %.4f s/km\n", amax, sheet.p * 1e3) * waveLabel;
        fontsize=13)
    return fig
end
devSave(recordSection(:x), "record_EW")
devSave(recordSection(:y), "record_NS")
nothing

# dev/japanTrench_nov2026/stages/step4_movie.jl
# Propagation movie (~20 s whatever the run length): left, the 3-D fence of u_z;
# right top, u_z on the rock top (coastline from GMT); right bottom, the E-W
# record section of u_z with a time cursor.

nFrames = length(snapTimes)
clim = robustLimit(3)
frame = Observable(1)
colours = @lift(fenceColours(3, $frame))

figMovie = Figure(size=(1700, 850))
axF = fenceAxis(figMovie[1:2, 1]; title=@lift(@sprintf("u_z, t = %5.1f s", snapTimes[$frame])))
hF = drawFence!(axF, @lift($colours[1]), @lift($colours[2]), @lift($colours[3]); clim=clim)
axM = Axis(figMovie[1, 2]; aspect=DataAspect(), xlabel="x, east (km)", ylabel="y, north (km)",
    title="u_z on the rock top (land surface / sea floor)")
heatmap!(axM, xKm, yKm, @lift($colours[1]); colormap=:balance, colorrange=(-clim, clim))
contour!(axM, xKm, yKm, bathymetry; levels=[0.0], color=:black, linewidth=1)
hlines!(axM, [yKm[jc]]; color=:gray, linestyle=:dash)
Colorbar(figMovie[1:2, 3], hF; label="u_z (m), clipped at the 99th percentile")

idxEW = 1:2:Nx
trEW = Float64.(surfaceU[3, idxEW, jc, :])
axR = Axis(figMovie[2, 2]; xlabel="time (s)", ylabel="x, east (km)", title="u_z along y = 0 (wiggles)")
wiggles!(axR, surfaceTimes, xKm[idxEW], trEW, 1.5 * 2dx / 1e3 / max(maximum(abs, trEW), eps()))
vlines!(axR, @lift([snapTimes[$frame]]); color=:red, linewidth=2)
xlims!(axR, 0, last(surfaceTimes))
colsize!(figMovie.layout, 1, Relative(0.52))
Label(figMovie[0, :], waveLabel; fontsize=14)

movieDir = joinpath(DEV_OUTPUT_DIR[], DEV_STEP[])
mkpath(movieDir)
moviePath = joinpath(movieDir, "planewave_$(japanTrenchCase).mp4")
framerate = max(1, round(Int, nFrames / 20))          # whole video ≈ 20 s
record(figMovie, moviePath, 1:nFrames; framerate=framerate) do n
    frame[] = n
end
println("movie -> ", moviePath, "  (", nFrames, " frames, ", framerate, " fps)")

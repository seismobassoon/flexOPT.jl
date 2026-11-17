# dev/japanTrench_nov2026/stages/step1_figures.jl
# Review figures of the box (outputs/step1_<case>/):
#   model_map        solid-top elevation, Vs of the top solid cell, NIED coverage
#   model_depth_*    Vp, Vs, NIED mask at DEPTH_SLICES_KM below sea level
#   model_section_EW Vp, Vs, ρ, Qμ along y = 0 (37°N);  model_section_NS along x = 0 (142°E)

xKm, yKm, zKm = xBox ./ 1e3, yBox ./ 1e3, zBox ./ 1e3
kSea = [argmin(abs.(@view altBox[i, j, :])) for i in 1:Nx, j in 1:Ny]
lonSurf = [geoBox[i, j, kSea[i, j]].lon for i in 1:Nx, j in 1:Ny]
latSurf = [geoBox[i, j, kSea[i, j]].lat for i in 1:Nx, j in 1:Ny]

"field with air -> NaN (and water -> NaN too if `water=false`)"
masked(field; water=true) = (f = Float64.(field); f[airMask] .= NaN; water || (f[waterMask] .= NaN); f)

function depthSlice(field, depthKm)
    out = fill(NaN, Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        k = argmin(abs.(@view(altBox[i, j, :]) .+ depthKm * 1e3))
        out[i, j] = field[i, j, k]
    end
    return out
end

"graticule (0.5°) + coastline + centre on a map axis"
function decorateMap!(ax)
    lonLevels = ceil(2minimum(lonSurf)) / 2:0.5:floor(2maximum(lonSurf)) / 2
    latLevels = ceil(2minimum(latSurf)) / 2:0.5:floor(2maximum(latSurf)) / 2
    contour!(ax, xKm, yKm, lonSurf; levels=collect(lonLevels), color=(:black, 0.35), linewidth=0.6)
    contour!(ax, xKm, yKm, latSurf; levels=collect(latLevels), color=(:black, 0.35), linewidth=0.6)
    contour!(ax, xKm, yKm, surfaceElevation; levels=[0.0], color=:black, linewidth=1.5)
    for L in lonLevels   # labels along the south edge
        i = findfirst(>=(L), @view lonSurf[:, 1])
        (isnothing(i) || i == 1) && continue
        t = (L - lonSurf[i-1, 1]) / (lonSurf[i, 1] - lonSurf[i-1, 1])
        text!(ax, xKm[i-1] + t * (xKm[i] - xKm[i-1]), yKm[1]; text=@sprintf("%.1f°E", L),
            align=(:center, :bottom), fontsize=9)
    end
    for L in latLevels   # labels along the west edge
        j = findfirst(>=(L), @view latSurf[1, :])
        (isnothing(j) || j == 1) && continue
        t = (L - latSurf[1, j-1]) / (latSurf[1, j] - latSurf[1, j-1])
        text!(ax, xKm[1], yKm[j-1] + t * (yKm[j] - yKm[j-1]); text=@sprintf("%.1f°N", L),
            align=(:left, :center), fontsize=9)
    end
    scatter!(ax, [0.0], [0.0]; marker=:star5, color=:red, markersize=14)
    return ax
end
mapAxis(pos, title) = Axis(pos; xlabel="x, east (km)", ylabel="y, north (km)", aspect=DataAspect(), title=title)

# ---------------------------------------------------------------- map
figMap = Figure(size=(1500, 520))
ax1 = mapAxis(figMap[1, 1], "solid-top elevation (land / sea floor)")
h1 = heatmap!(ax1, xKm, yKm, surfaceElevation ./ 1e3; colormap=:oleron, colorrange=(-8, 8))
decorateMap!(ax1); Colorbar(figMap[2, 1], h1; vertical=false, label="km")
ax2 = mapAxis(figMap[1, 2], "Vs of the top solid cell")
topVs = [VsModel[i, j, topRock[i, j]] for i in 1:Nx, j in 1:Ny]
h2 = heatmap!(ax2, xKm, yKm, topVs; colormap=:viridis)
decorateMap!(ax2); Colorbar(figMap[2, 2], h2; vertical=false, label="km/s")
ax3 = mapAxis(figMap[1, 3], "NIED fraction of the rock column (rest = PREM)")
h3 = heatmap!(ax3, xKm, yKm, niedColumnFraction; colormap=:grays, colorrange=(0, 1))
decorateMap!(ax3); Colorbar(figMap[2, 3], h3; vertical=false, label="fraction")
Label(figMap[0, :], @sprintf("Japan Trench box, centre %.1f°N %.1f°E, dx = %.1f km (%s)",
    BOX_CENTRE.lat, BOX_CENTRE.lon, dx / 1e3, japanTrenchCase); fontsize=16)
devSave(figMap, "model_map")

# ---------------------------------------------------------------- depth slices
VpPlot, VsPlot = masked(VpModel; water=false), masked(VsModel)
vpRange = extrema(VpModel[rockMask]); vsRange = (VS_MIN_SOLID, maximum(VsModel[rockMask]))
nD = length(DEPTH_SLICES_KM)
figD = Figure(size=(320nD + 80, 1000))
for (c, d) in enumerate(DEPTH_SLICES_KM)
    a = mapAxis(figD[1, c], @sprintf("Vp, %g km", d))
    hp = heatmap!(a, xKm, yKm, depthSlice(VpPlot, d); colormap=:turbo, colorrange=vpRange)
    decorateMap!(a)
    b = mapAxis(figD[2, c], @sprintf("Vs, %g km (blue = water)", d))
    hs = heatmap!(b, xKm, yKm, depthSlice(VsPlot, d); colormap=:turbo, colorrange=vsRange, lowclip=:royalblue)
    decorateMap!(b)
    m = mapAxis(figD[3, c], @sprintf("NIED used, %g km", d))
    hm = heatmap!(m, xKm, yKm, depthSlice(Float64.(niedMask .& rockMask), d); colormap=:grays, colorrange=(0, 1))
    decorateMap!(m)
    if c == nD
        Colorbar(figD[1, nD+1], hp; label="Vp (km/s)")
        Colorbar(figD[2, nD+1], hs; label="Vs (km/s)")
        Colorbar(figD[3, nD+1], hm; label="1 = NIED, 0 = PREM")
    end
end
Label(figD[0, :], "depth slices below sea level (" * japanTrenchCase * ")"; fontsize=16)
devSave(figD, "model_depth_slices")

# ---------------------------------------------------------------- sections
ρPlot, QμPlot = masked(ρModel; water=false), masked(QμModel; water=false)
function sectionFigure(along, name)
    sec(f) = along == :x ? f[:, jc, :] : f[ic, :, :]
    h = along == :x ? xKm : yKm
    lab = along == :x ? "x, east (km)" : "y, north (km)"
    fig = Figure(size=(1300, 1200))
    rows = [("Vp (km/s)", sec(VpPlot), :turbo, vpRange, nothing),
            ("Vs (km/s), blue = water", sec(VsPlot), :turbo, vsRange, :royalblue),
            ("ρ (g/cm³)", sec(ρPlot), :viridis, extrema(ρModel[rockMask]), nothing),
            ("Qμ (planet1D)", sec(QμPlot), :plasma, extrema(QμModel[rockMask]), nothing)]
    for (r, (t, f, cm, cr, lc)) in enumerate(rows)
        ax = Axis(fig[r, 1]; xlabel=lab, ylabel="z (km, local Cartesian)", title=t, aspect=4.0)
        hm = isnothing(lc) ? heatmap!(ax, h, zKm, f; colormap=cm, colorrange=cr) :
             heatmap!(ax, h, zKm, f; colormap=cm, colorrange=cr, lowclip=lc)
        Colorbar(fig[r, 2], hm)
    end
    Label(fig[0, :], name * @sprintf("  (sea level is %.1f km below z = 0 at the box edge)",
        (earthRadius - sqrt(earthRadius^2 - BOX_HALF_WIDTH^2)) / 1e3); fontsize=16)
    return fig
end
devSave(sectionFigure(:x, @sprintf("E-W section through %.1f°N", BOX_CENTRE.lat)), "model_section_EW")
devSave(sectionFigure(:y, @sprintf("N-S section through %.1f°E", BOX_CENTRE.lon)), "model_section_NS")
nothing

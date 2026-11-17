# dev/japanTrench_nov2026/stages/step1_box.jl
# Local Cartesian box (GeoPoints.constructLocalBox). The box is a true
# Cartesian frame: at the corners (283 km from the centre) the sea surface is
# ~6 km below the tangent plane, so the bottom is lowered by that much to keep
# DEPTH_BELOW_SEA_LEVEL of Earth under every column.

p0 = GeoPoint(BOX_CENTRE.lat, BOX_CENTRE.lon)
earthRadius = planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer * 1e3
cornerDistance = sqrt(2) * BOX_HALF_WIDTH
cornerDrop = earthRadius - sqrt(earthRadius^2 - cornerDistance^2)
altMin = -ceil((DEPTH_BELOW_SEA_LEVEL + cornerDrop) / dx) * dx
altMax = ceil((TOPO_MAX + dx) / dx) * dx          # highest peak + >= 1 air layer

box = constructLocalBox(
    p0, dx, dx, dx,
    -BOX_HALF_WIDTH, BOX_HALF_WIDTH,
    -BOX_HALF_WIDTH, BOX_HALF_WIDTH,
    altMin, altMax;
    axis_angle_deg=BOX_AXIS_AZIMUTH,
)
geoBox = box.allGridsInGeoPoints
Nx, Ny, Nz = size(geoBox)
xBox = [p.xyz[1] for p in box.allGridsInCartesian[:, 1, 1]]
yBox = [p.xyz[2] for p in box.allGridsInCartesian[1, :, 1]]
zBox = [p.xyz[3] for p in box.allGridsInCartesian[1, 1, :]]
altBox = [p.alt for p in geoBox]                  # m above the WGS84 ellipsoid (≈ sea level)
@assert isapprox(xBox[2] - xBox[1], dx) && isapprox(zBox[2] - zBox[1], dx)

# orientation and coverage checks
ic, jc = argmin(abs.(xBox)), argmin(abs.(yBox))
@assert geoBox[end, jc, Nz].lon > geoBox[1, jc, Nz].lon "x should point east"
@assert geoBox[ic, end, Nz].lat > geoBox[ic, 1, Nz].lat "y should point north"
@assert maximum(@view altBox[:, :, 1]) <= -DEPTH_BELOW_SEA_LEVEL + 1.0 "bottom too shallow somewhere"
@assert minimum(@view altBox[:, :, end]) >= TOPO_MAX "top too low somewhere"

@printf("box: %d x %d x %d points (%.2f M), dx = %.1f km\n", Nx, Ny, Nz, Nx * Ny * Nz / 1e6, dx / 1e3)
@printf("     x (east), y (north) in [%.0f, %.0f] km; z in [%.0f, %.0f] km; sea level %.1f km below the tangent plane at the corners\n",
    first(xBox) / 1e3, last(xBox) / 1e3, first(zBox) / 1e3, last(zBox) / 1e3, cornerDrop / 1e3)
for (name, (i, j)) in (("SW", (1, 1)), ("SE", (Nx, 1)), ("NE", (Nx, Ny)), ("NW", (1, Ny)))
    p = geoBox[i, j, Nz]
    @printf("     %s corner: %.3f°N %.3f°E\n", name, p.lat, p.lon)
end

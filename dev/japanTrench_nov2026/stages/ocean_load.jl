# dev/japanTrench_nov2026/stages/ocean_load.jl
# Ocean as a load ("burden") on the sea floor.
#
# NB getParamsAndTopo (src/GeoPoints/getSeismicParamTopo.jl) tests "above the
# topography -> air" BEFORE "between sea floor and sea level -> water", so the
# ocean comes out as AIR: the top rock cell of an oceanic column is already the
# sea floor, i.e. a free surface there. Here we only add the water's mass:
#   water depth h(lon, lat) from the same GMT earth_relief grid (cache hit in
#   data/topoViaGMT), ρ_top += ρw h / dx on the top rock cell, λ and μ kept.
# The added mass acts on all three components (the OPT operator has one ρ per
# node); physically an inviscid ocean loads the vertical motion only, so the
# horizontal inertia of the oceanic sea floor is overestimated.

topoRegion = [extrema(p.lon for p in geoBox)..., extrema(p.lat for p in geoBox)...]
topoGrid = let precision = GeoPoints.GMTprecision(TOPO_PRECISION_KM), region = topoRegion
    myProduceOrLoad(GeoPoints.getTopoViaGMT, @strdict(precision, region), "topoViaGMT")["topo"]
end
topoItp = let gx = topoGrid.x, gy = topoGrid.y, gz = topoGrid.z
    if size(gx, 1) != size(gz, 2)       # same pixel/grid-registration fix as getParamsAndTopo
        gx = (gx[1:end-1] .+ gx[2:end]) ./ 2
        gy = (gy[1:end-1] .+ gy[2:end]) ./ 2
    end
    linear_interpolation((gx, gy), permutedims(gz); extrapolation_bc=Flat())
end

kSea = [argmin(abs.(@view altBox[i, j, :])) for i in 1:Nx, j in 1:Ny]
lonSurf = [geoBox[i, j, kSea[i, j]].lon for i in 1:Nx, j in 1:Ny]
latSurf = [geoBox[i, j, kSea[i, j]].lat for i in 1:Nx, j in 1:Ny]
bathymetry = [topoItp(lonSurf[i, j], latSurf[i, j]) for i in 1:Nx, j in 1:Ny]   # m, < 0 at sea
waterDepth = max.(0.0, .-bathymetry)

# the free surface is the top of the rock: water (if a later getParamsAndTopo
# ever keeps it) is treated as air here too
materialFS = copy(rockMask)
ρLoaded = copy(ρModel)                    # g/cm³; only ρ of the top rock cells changes
if OCEAN_LOAD === :top_density
    for j in 1:Ny, i in 1:Nx
        waterDepth[i, j] > 0 || continue
        ρLoaded[i, j, topRock[i, j]] += WATER.ρ * waterDepth[i, j] / dx
    end
elseif OCEAN_LOAD !== :none
    error("OCEAN_LOAD must be :top_density or :none")
end
@printf("ocean load (%s): %.0f %% of the columns are at sea, water depth up to %.2f km; top-cell ρ up to %.2f g/cm³ (rock alone ≤ %.2f)\n",
    OCEAN_LOAD, 100count(waterDepth .> 0) / length(waterDepth), maximum(waterDepth) / 1e3,
    maximum(ρLoaded[i, j, topRock[i, j]] for i in 1:Nx, j in 1:Ny),
    maximum(ρModel[i, j, topRock[i, j]] for i in 1:Nx, j in 1:Ny))

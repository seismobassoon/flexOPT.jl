# dev/japanTrench_nov2026/stages/config.jl
#
# Every choice of the Japan Trench box, in one place. `japanTrenchCase`
# (quick | standard | full) only changes the grid spacing.

# ---------------------------------------------------------------- box
# Centre 37°N 142°E (off Fukushima / Miyagi), ±200 km east and north.
# GeoPoints.constructLocalBox: axis_angle_deg = azimuth of x (0 = N, 90 = E);
# y = z × x, so 90° gives x -> east, y -> north, z -> up.
BOX_CENTRE = (lat=37.0, lon=142.0)
BOX_AXIS_AZIMUTH = 90.0
BOX_HALF_WIDTH = 200e3            # m, in x (east) and y (north)
DEPTH_BELOW_SEA_LEVEL = 80e3      # m: every column, corners included, reaches 80 km below sea level
TOPO_MAX = 3.0e3                  # m kept above sea level (highest peaks here ~2 km: Azuma, Bandai, Nasu)

# ---------------------------------------------------------------- grids
# Isotropic grid (dx = dy = dz). minPeriod ≈ 4 dx / VS_MIN_SOLID (OPT3, ~4
# points per S wavelength in the slowest solid). The paper works at 20-50 s.
GRID_PRESETS = Dict(
    "coarse"   => (dx=8.0e3, minPeriod=32.0),
    "quick"    => (dx=5.0e3, minPeriod=20.0),
    "standard" => (dx=2.5e3, minPeriod=10.0),
    "full"     => (dx=1.0e3, minPeriod=4.0),
)
preset = GRID_PRESETS[japanTrenchCase]
dx = preset.dx
minResolvedPeriod = preset.minPeriod

# ---------------------------------------------------------------- materials
# NIED 3-D model (same archive and same settings as the Kirishima box):
# low-confidence or outside points fall back to planet1D (PREM), which also
# gives ρ, Qμ and Qκ everywhere (NIED has Vp and Vs only).
niedVelocitySource = DEFAULT_NIED_VELOCITY_SOURCE[]
velocitySourceKey = nied_velocity_source_key(niedVelocitySource)
NIED_CONFIDENCE_MAX = 0.8
NIED_OUTSIDE = :planet1D
NIED_LOW_CONFIDENCE = :planet1D
TOPO_PRECISION_KM = 1.0           # GMTprecision -> @earth_relief_15s
WATER = (ρ=1.03, vp=1.5)          # g/cm³, km/s; the paper: constant Vp = 1.5 km/s sea water
AIR = (ρ=0.001, vp=0.314)         # getParamsAndTopo defaults
AIR_DENSITY_CUTOFF = 0.01         # g/cm³, rock/water vs air (Kirishima value)
VS_MIN_SOLID = 1.0                # km/s, as Yamaya et al. (2024) §3 (applied after the cache, Vp kept)

# ---------------------------------------------------------------- step2: OPT3 propagation
# Cerjan layers on the four sides and the bottom (none on top: free surface).
CERJAN_CELLS = Dict("coarse" => 10, "quick" => 12, "standard" => 16, "full" => 20)
cerjanCells = CERJAN_CELLS[japanTrenchCase]
cerjanDamping = 3.0 / cerjanCells^2   # edge factor exp(-3), Kirishima/Greenland value

# Ocean "as a burden": free surface on the sea floor (getParamsAndTopo already
# makes the water air, see README) and the water column's mass per unit area,
# ρw h, added to the top rock cell: ρ_top += ρw h / dx (λ, μ unchanged).
# :none = no load (sea floor is a plain free surface).
OCEAN_LOAD = :top_density

# Plane wave from a smooth horizontal "sheet" at PLANE_WAVE.depth (Cartesian
# z = -depth): moment density m(x,t) = M w(x,y) G(z) r(t - p·x_h), with
#   P : f = -∇m      (curl-free  -> P only, in a homogeneous layer)
#   SV: f = ∇m × e⊥  (div-free   -> S, polarised in the propagation plane)
#   SH: f = ∇m × ẑ   (div-free   -> S, horizontal)
# p = sin(incidence)/V(depth) (median Vp or Vs of that layer), r = Ricker of
# peak period `period`, w = Tukey taper over `taper` m inside the box edges,
# G = Gaussian of half-width 1.5 dx. Both an up- and a down-going wave leave
# the sheet; the down-going one goes into the bottom Cerjan layer.
PLANE_WAVE = (
    type=:P,              # :P | :SV | :SH
    incidence=45.0,       # deg from the vertical, at the sheet depth
    azimuth=270.0,        # deg from N, direction of horizontal propagation (270: towards W, i.e. from the trench side)
    period=20.0,          # s, Ricker peak period
    depth=65e3,           # m below z = 0 (Cartesian)
    taper=60e3,           # m, lateral cosine taper inside the physical box
    amplitude=1e6,        # Pa, moment density scale (everything is linear in it)
)
duration = 160.0          # s
courant = 0.25            # dt = courant * dx / (sqrt(3) * Vp_max)
optTimeMarchingScheme = :predictor_corrector
optPredictorMode = :lumped_future
optCorrectorIterations = 8
optCorrectorTolerance = 1e-6
optCorrectorRelaxation = 1.0
blowupLimit = 1e3         # m
traceSampling = 0.5       # s, surface wavefield (every rock-top node, 3 components)
snapshotEvery = 2.0       # s, vertical E-W and N-S sections (multiple of traceSampling)
optSupplementaryOrder = 2 # OPT3, as Kirishima / Greenland (same recipe cache)
recipeOnGPU = false

# ---------------------------------------------------------------- figures
DEPTH_SLICES_KM = [5.0, 15.0, 30.0, 50.0, 70.0]   # below sea level

JAPANTRENCH_CACHE_DIR = joinpath("dev", "japanTrench_nov2026")
println("Japan Trench case = ", japanTrenchCase, ": dx = ", dx / 1e3,
    " km, periods ≳ ", minResolvedPeriod, " s resolved with Vs ≥ ", VS_MIN_SOLID, " km/s")

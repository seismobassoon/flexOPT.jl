# dev/new_Paul_models/stages/config.jl -- edit freely, re-read on every run
PAUL_GRID_N = 701            # Cartesian pixels across the disk
PAUL_SEAM_BLEND_DEG = 0.0    # 0 = raw StagYY data (batchStagYY default elsewhere is 2°)
PAUL_SKIP_CODES = String[]   # e.g. ["phs11"] to skip some files

# Pseudo nₑ (figure <snap>_01_electron_density): end-member compositions of
# Paul's runs (mol%, NCFMAS, from Paul 2026-11-15) -> Yₑ = Σ xᵢZᵢ / Σ xᵢAᵢ
# (harzburgite 0.495573, basalt 0.495049). bs is a MASS fraction (Paul).
PAUL_COMP_HZ = Dict("SiO2" => 36.8777, "MgO" => 55.9288, "FeO" => 5.9150,
    "CaO" => 0.5749, "Al2O3" => 0.7017, "Na2O" => 0.0019)
PAUL_COMP_MORB = Dict("SiO2" => 47.5554, "MgO" => 20.5098, "FeO" => 7.3961,
    "CaO" => 14.3570, "Al2O3" => 9.5499, "Na2O" => 0.6318)
PAUL_YE_HZ = yeFromMolPercent(PAUL_COMP_HZ)
PAUL_YE_MORB = yeFromMolPercent(PAUL_COMP_MORB)
PAUL_BS_IS_VOLUME = false     # bs is a mass fraction
PAUL_DRHO_MORB = 0.02         # only used if PAUL_BS_IS_VOLUME = true
@printf("Yₑ harzburgite = %.6f, basalt = %.6f\n", PAUL_YE_HZ, PAUL_YE_MORB)

# ---- synthetics (step `synthetics`) and phase windows (step `phases`) --------
# FD3 whole-Earth 2-D P-SV disk (as scripts/earth2D_StagYY_neutrino_seismo.jl).
SEISMO = (
    dx = 10e3,                 # m
    period = 30.0,             # s, moment-rate spectrum exp(−(f T)²)
    duration = 2400.0,         # s
    cfl = 0.38,
    surfaceDamping = 2e-3,     # per-step damping on the outermost material cells only
    outputSampling = 1.0,      # s
    snapshotTimes = (300.0, 600.0, 900.0, 1200.0),
    mohoRadius = 6346.6e3,     # LID above (no crust, no ocean)
    clip = 0.10,               # |δlnV|, |δlnρ| ≤ clip
    sourceDepth = 2000e3,      # m, both sources (north and south pole)
    M0 = 1e15,                 # N m per m of out-of-plane length (linear: only a scale)
    dip = 45.0, slip = :thrust,# in-plane dip-slip: Mrr = M0, Mtt = −M0, Mrt = 0 for 45° thrust
    receiverStepDeg = 2.0,     # stations every 2° all around (180 stations)
    receiverDepth = 20e3,      # = 2 dx, below the staircase surface
)
SEISMO_MODELS = [(model = "G1m20_sb11", snap = "03140"), (model = "G1m20_sb21", snap = "03240")]
SEISMO_POLES = (:north, :south)
SEISMO_CACHE_DIR = joinpath("dev", "new_Paul_models")     # data/dev/new_Paul_models/*.jld2

# Phases (names of dev/raytheory/lib/rayPhases.jl). Last leg P -> vertical
# component, last leg S -> horizontal (radial, SV: no SH in 2-D P-SV).
SEISMO_PHASES = ["P", "pP", "PcP", "PP", "PKP", "PKIKP", "S", "sS", "ScS", "SKS", "SS"]
SEISMO_WINDOW = (pre = 25.0, post = 45.0, maxLag = 15.0)  # s around the PREM time
SEISMO_SECTION_PHASES = ["P", "pP", "PcP", "PP", "PKP", "S", "sS", "ScS", "SKS", "SS"]

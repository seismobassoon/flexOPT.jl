# dev/new_Paul_models/stages/neutrino_env.jl -- neutrino modules (as scripts/earth2D_StagYY_neutrino_seismo.jl)
nuLoadEnv = !@isdefined(PAUL_NU_ENV)
nuLoadEnv && !@isdefined(neutrinoOscillation) && include(joinpath(flexopt_root, "src", "neutrinoOscillation.jl"))
nuLoadEnv && !@isdefined(GeoPoints) && include(joinpath(flexopt_root, "src", "GeoPoints.jl"))
nuLoadEnv && !@isdefined(neutrinoFlux) && include(joinpath(flexopt_root, "src", "neutrinoFlux.jl"))
nuLoadEnv && !@isdefined(Neutrino_Cross_Sections) && include(joinpath(flexopt_root, "src", "neutrinoCrossSections.jl"))
nuLoadEnv && !@isdefined(neutrinoDetector) && include(joinpath(flexopt_root, "src", "neutrinoDetector.jl"))
using .neutrinoOscillation, .GeoPoints, .neutrinoFlux, .Neutrino_Cross_Sections, .neutrinoDetector
using StaticArrays
PAUL_NU_ENV = true

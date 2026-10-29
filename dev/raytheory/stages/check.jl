# dev/raytheory/stages/check.jl — every numerical check, no plotting package needed.
# Cases (models for the spherical part): default PREM; e.g. `check PREM ak135 DWThot vprem`.
allChecks = RTCheck[]
append!(allChecks, printChecks(checkCartesian(); title="Cartesian, isotropic"))
for case in (isempty(DEV.cases) ? ["PREM"] : DEV.cases)
    m = rtModel(case)
    show(stdout, m)
    append!(allChecks, printChecks(checkSphericalTracer(m); title="$case: tracer vs τ-p (same ray parameter)"))
    case == "PREM" && append!(allChecks, printChecks(checkTauP(m, RT_TAUP_CSV); title="PREM: τ-p vs TauP"))
end
append!(allChecks, printChecks(checkAnisotropy(); title="anisotropy (Christoffel Hamiltonian)"))
@printf("\nALL: %d / %d checks passed\n", count(passed, allChecks), length(allChecks))

# Outputs_1D — derived quantities with OPT precision ("cooking recipes")

The solver gives a field u, but what we want is often a quantity *derived*
from it by a PDE expression: strain, stress, flux, traction, divergence…
Here such a quantity is **declared exactly like a famous equation** and
turned into an OPT recipe by the same machinery (`makeOPTsemiSymbolic`) and
the same scheme parameters as the operator:

| entry | expression | role |
|---|---|---|
| `1DoutFlux` | q = a(x) ∂x u | the derived quantity (1-D analogue of σ = μ ∂x u) |
| `1DoutValue` | m = c(x) u, c = 1 | the test-function "mass", to turn weak-form rows into point values |

Test: u = cos x, a(x) = 2 + sin(x/3), exact q = −a sin x. Compared methods:

* **FD3 / FD5 (Fornberg)**: what `dev/spherical1DelasticKernel` uses today for strains;
* **⟨scheme⟩ lumped**: q_ν = (B u)_ν / (M 1)_ν with B, M the recipe matrices;
* **⟨scheme⟩ consistent**: q = M⁻¹ B u (Petrov–Galerkin projection).

for the schemes convFD3, OPT3, OPT5. Boundary rows are replaced by one-sided
FD5 rows and excluded from the error (12 nodes at each end, shaded in the
figure); the free-surface / traction rows are the next step.

## Run (on stella: inside `srun` or a job, never on `stella*`)

```sh
julia --project=. application/Outputs_1D/run.jl                     # compare, G = 8 points per wavelength
julia --project=. application/Outputs_1D/run.jl compare G=5
julia --project=. application/Outputs_1D/run.jl convergence         # error vs G = 4 … 40
julia --project=. application/Outputs_1D/run.jl convergence schemes=OPT3,OPT5
```

Figures and CSV: `application/Outputs_1D/outputs/compare/`, `outputs/convergence/`.
The 6 recipes (2 entries × 3 schemes, built at Δ = 1) are cached in `data/semiSymbolic/`.

## Where this goes

1. this application: 1-D Cartesian flux (step 1);
2. 1-D spherical SH output σ_rφ;
3. P-SV ε, σ with free-surface and traction-continuity rows (homogeneous sphere, Lamb);
4. then a `famousOutputs` gallery in `src/motorsOPT/` and `computeSGTDisk` using the output rows
   for strains and their transpose for moment-tensor sources (exact reciprocity kept).

Nothing in `src/` is modified: the entries are new methods of `flexOPT.famousEquation`.

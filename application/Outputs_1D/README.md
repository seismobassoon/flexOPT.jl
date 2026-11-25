# Outputs_1D — derived quantities chosen as PDE expressions, cooked with OPT recipes

The solver gives a field u, but what we want is often a quantity *derived* from it: strain,
stress, flux, traction, divergence… Here **the user writes that quantity as a PDE expression**
and flexOPT turns it into an OPT recipe with the same machinery (`makeOPTsemiSymbolic`) and the
same scheme parameters as the operators:

```julia
registerOutput!("stress", (u, μ) -> μ * ∂x(u))                    # q = μ ∂x u
registerOutput!("mixed",  (u, a, b) -> a * ∂x(u) + b * u)          # any terms, any orders
lumped, consistent = cookOutput("stress", uvals, [μvals], Δx; scheme="OPT3", mass="Gal3")
```

`registerOutput!` creates the famousEquation entry `out_<name>` (the expression is part of the
recipe cache key). Point values come from a projection, q = M⁻¹ B u: B is the output recipe and
M the recipe of the "value" output c·u (c = 1), from the same scheme or another one (`mass=`).
"lumped" is (B u)/(M 1), without the solve.

## First result (output a ∂x u, G = 8 points per wavelength)

| method | error |
|---|---|
| convFD3 (= FD3) | 10 % |
| convFD5 (= FD5) | 1.2 % |
| OPT3 consistent | 5.4 % |
| OPT5 consistent | 0.59 % |
| convFD3 consistent | 0.31 % (the classical 4th-order Padé compact derivative) |

OPT schemes modify their mass to cancel the errors of the wave operator ω²M − K; for an output
projection that is the wrong partner. This version therefore also tries the OPT output rows B
with a plain Petrov–Galerkin mass (`Gal3`, `Gal5` = OPT3/OPT5 parameters with
`supplementaryOrder = 0`, experimental).

## Run (on stella: inside `srun` or a job, never on `stella*`)

```sh
julia --project=. application/Outputs_1D/run.jl                              # compare, output "flux", G = 8
julia --project=. application/Outputs_1D/run.jl compare name=mixed G=6
julia --project=. application/Outputs_1D/run.jl convergence                  # G = 4 … 32
julia --project=. application/Outputs_1D/run.jl custom expr="a*∂x(∂x(u))+b*∂x(u)"
julia --project=. -e 'include("application/Outputs_1D/examples.jl")'   # four worked examples
```

`custom`: field `u`, materials `a`, `b`, `c` (= 2 + sin(x/3), 1 + cos(x/5)/2, 1.5; at least one
must appear), no commas in the expression. Test field u = cos x on [0, 60]; the exact q is
computed symbolically from the same expression. Pairs compared (B scheme, M scheme): convFD3,
convFD5, OPT3, OPT5, OPT3 + Gal3 mass, OPT5 + Gal5 mass, Gal3. A pair whose recipe fails is
skipped with a warning.

Recipes are built at the real spacing Δ (needed for expressions mixing derivative orders), so a
convergence run builds a few recipes per G the first time (then cached in `data/semiSymbolic/`).
Boundary rows: the first/last `pointsInSpace` nodes are set to the exact value and excluded
from the error (12 nodes at each end); free-surface / traction rows are the next step.

## Where this goes

1. this application: 1-D Cartesian outputs, user expressions (now);
2. 1-D spherical SH output σ_rφ;
3. P-SV ε, σ with free-surface and traction-continuity rows (homogeneous sphere, Lamb);
4. a `famousOutputs` gallery in `src/motorsOPT/`, and `computeSGTDisk` using output rows for the
   strains and their transpose for moment-tensor sources (exact reciprocity kept).

Nothing in `src/` is modified: the entries are new methods of `flexOPT.famousEquation`.

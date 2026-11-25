# application/Outputs_1D/examples.jl -- choose your own output expressions (REPL, one statement per line)
#
#   julia --project=.
#   julia> include("application/Outputs_1D/examples.jl")
# or line by line after include("application/Outputs_1D/run.jl").
# Each registerOutput! must be its own top-level statement before the output is used (world age).

isdefined(Main, :registerOutput!) || include(joinpath(@__DIR__, "run.jl"))

# 1. stress-like flux: q = μ ∂x u
registerOutput!("stress", (u, μ) -> μ * ∂x(u))
r1 = outputs1DCompare("stress"; G=6.0)

# 2. several terms, mixed derivative orders: q = a ∂x u + b u
registerOutput!("mixed2", (u, a, b) -> a * ∂x(u) + b * u)
r2 = outputs1DCompare("mixed2"; G=6.0)

# 3. a second-order output (e.g. a curvature / bending term): q = a ∂x² u + b ∂x u
registerOutput!("curvature", (u, a, b) -> a * ∂x(∂x(u)) + b * ∂x(u))
r3 = outputs1DCompare("curvature"; G=8.0)

# 4. your own data: cook an output from nodal values (here a Gaussian bump) with OPT3 rows + Galerkin mass
Δx = 0.1
X = collect(0.0:Δx:20.0)
uvals = exp.(-(X .- 10.0) .^ 2)
μvals = 2.0 .+ sin.(X ./ 3)
lumped, consistent = cookOutput("stress", uvals, [μvals], Δx; scheme="OPT3", mass="Gal3")
exact = μvals .* (-2 .* (X .- 10.0) .* uvals)
println("example 4 (Gaussian): max |consistent - exact| = ", maximum(abs.(consistent[10:end-9] .- exact[10:end-9])))

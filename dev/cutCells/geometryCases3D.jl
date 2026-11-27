# dev/cutCells/geometryCases3D.jl -- the 3-D topography cases, defined ONCE for every code
#
# Plain Julia (LinearAlgebra only). Coordinates (x, y, z), z positive upward, metres.
# Same medium, source time function and grid spacing as the 2-D cases (geometryCases.jl);
# the box is smaller (3-D cost).
#
# Cases
#   flat3D          free surface z = 0                                  (θ = 0)
#   oblique10       plane with unit normal n = (-sinθ cosφ, -sinθ sinφ, cosθ), θ = 10°, φ = 30°
#   oblique20       same, θ = 20°, φ = 30°
# An oblique case is the flat case RIGIDLY ROTATED by R (rotation about ez × n by θ):
# surface, source, force (= n) and receivers rotate together, the grid does not.
# So its components (u_1, u_2, u_n) = Rᵀ u equal the flat (u_x, u_y, u_z) at the same
# receiver: the exact reference is the flat solution.
#
# Flat-frame geometry: source (0, 0, -6 km), force along +z (outward normal);
# receivers on two lines, 9 each: x' = -5, -2.5, 0, 2.5, 5 km at y' = 0, then
# y' = -5, -2.5, 2.5, 5 km at x' = 0; deep line z' = -3 km, shallow line z' = -1 km.
# Box: ±8 km around the source, 8 km below it, 10 Cerjan cells (sides, bottom), 6 s --
# sized for memory (OPT3 3-D operators ~15 kB per node, plus copies while installing rows).
# Receiver x'=y'=0 is directly above the source: u_1 = u_2 = 0 there by symmetry.

using LinearAlgebra

const GC3 = (vp=6000.0, vs=3464.0, rho=2700.0, dx=500.0,
    duration=6.0, force=1e10, sigma=500.0, cutoff=2.0,
    sourceDepth=6e3, deepDepth=3e3, shallowDepth=1e3,
    lineOffsets=[-5e3, -2.5e3, 0.0, 2.5e3, 5e3], crossOffsets=[-5e3, -2.5e3, 2.5e3, 5e3],
    halfWidth=8e3, below=8e3, voidCapCells=3, cerjanCells=10, azimuth=30.0)

gc3Frequency() = 1.5 * GC3.vs / (16 * 500.0)          # 0.649500 Hz, as in 2-D
gc3Delay() = 1.5 / gc3Frequency()                    # 2.309469 s
function gc3Ricker(t)
    a = π * gc3Frequency() * (t - gc3Delay())
    (1 - 2a^2) * exp(-a^2)
end

const CASES3D = ("flat3D", "oblique10", "oblique20")
case3DTheta(case) = case == "flat3D" ? 0.0 : case == "oblique10" ? 10.0 :
    case == "oblique20" ? 20.0 : error("unknown 3-D case $case")
case3DAzimuth(case) = GC3.azimuth

"Unit outward normal of the free surface."
function case3DNormal(case)
    θ, φ = deg2rad(case3DTheta(case)), deg2rad(case3DAzimuth(case))
    [-sin(θ) * cos(φ), -sin(θ) * sin(φ), cos(θ)]
end

"Rotation flat frame -> grid frame (R ez = n); identity for flat3D."
function case3DRotation(case)
    n = case3DNormal(case)
    ez = [0.0, 0.0, 1.0]
    a = cross(ez, n)
    s = norm(a)
    s < 1e-14 && return Matrix{Float64}(I, 3, 3)
    a ./= s
    c = dot(ez, n)
    K = [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
    Matrix{Float64}(I, 3, 3) + s * K + (1 - c) * K * K
end

case3DToGrid(case, p) = case3DRotation(case) * collect(Float64, p)
"Grid components -> (u_1, u_2, u_n) (flat-frame components)."
case3DRotate(case, ux, uy, uz) = begin
    Rt = permutedims(case3DRotation(case))
    (Rt[1, 1] .* ux .+ Rt[1, 2] .* uy .+ Rt[1, 3] .* uz,
     Rt[2, 1] .* ux .+ Rt[2, 2] .* uy .+ Rt[2, 3] .* uz,
     Rt[3, 1] .* ux .+ Rt[3, 2] .* uy .+ Rt[3, 3] .* uz)
end

"Level set φ(x, y, z) = n · r: material where φ ≤ 0."
case3DPhi(case) = (x, y, z) -> dot(case3DNormal(case), (x, y, z))
"Surface height z_s(x, y)."
function case3DSurface(case)
    n = case3DNormal(case)
    (x, y) -> -(n[1] * x + n[2] * y) / n[3]
end
case3DSource(case) = case3DToGrid(case, (0.0, 0.0, -GC3.sourceDepth))
case3DForce(case) = case3DNormal(case)

"Receivers: Dict(:deep|:shallow => [(x', y', grid point)]) in the order of the header."
function case3DReceivers(case)
    flat = vcat([(o, 0.0) for o in GC3.lineOffsets], [(0.0, o) for o in GC3.crossOffsets])
    Dict(line => [(xp=p[1], yp=p[2], r=case3DToGrid(case, (p[1], p[2], -depth))) for p in flat]
         for (line, depth) in ((:deep, GC3.deepDepth), (:shallow, GC3.shallowDepth)))
end

"""
    case3DGrid(case; dx=GC3.dx, padding=GC3.cerjanCells, halfWidth=GC3.halfWidth, below=GC3.below)

Model box (x, y, z axes, node multiples of dx) centred laterally on the source,
GC3.halfWidth to each side and GC3.below under the source, up to the highest surface
point over the PADDED lateral extent plus GC3.voidCapCells. Returns
(x, y, z) of the model box and (xw, yw, zw) of the whole grid including the
Cerjan padding (sides and bottom, none on top).
"""
function case3DGrid(case; dx=GC3.dx, padding=GC3.cerjanCells, halfWidth=GC3.halfWidth, below=GC3.below)
    s = case3DSource(case)
    snap(v) = round(v / dx) * dx
    cx, cy, cz = snap(s[1]), snap(s[2]), snap(s[3])
    W = halfWidth
    x = collect((cx - W):dx:(cx + W))
    y = collect((cy - W):dx:(cy + W))
    P = padding * dx
    zs = case3DSurface(case)
    top = maximum(zs(xx, yy) for xx in (x[1] - P, x[end] + P), yy in (y[1] - P, y[end] + P))
    z = collect((cz - below):dx:(ceil(top / dx) * dx + GC3.voidCapCells * dx))
    xw = collect((x[1] - P):dx:(x[end] + P))
    yw = collect((y[1] - P):dx:(y[end] + P))
    zw = collect((z[1] - P):dx:z[end])
    (x=x, y=y, z=z, xw=xw, yw=yw, zw=zw)
end

"Gaussian source weights on the nodes of (x, y, z): [(i, j, k, w)], Σw = 1."
function case3DSourceWeights(case, x, y, z)
    s = case3DSource(case)
    σ, rc = GC3.sigma, GC3.cutoff * GC3.sigma
    out = Tuple{Int,Int,Int,Float64}[]
    for k in eachindex(z), j in eachindex(y), i in eachindex(x)
        d2 = (x[i] - s[1])^2 + (y[j] - s[2])^2 + (z[k] - s[3])^2
        d2 <= rc^2 && push!(out, (i, j, k, exp(-d2 / (2σ^2))))
    end
    total = sum(last, out)
    [(i, j, k, w / total) for (i, j, k, w) in out]
end

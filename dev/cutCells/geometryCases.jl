# dev/cutCells/geometryCases.jl -- ONE definition of every topography benchmark case,
# shared by flexOPT (tiltedPlane.jl), SPECFEM2D (specfemCases.jl) and SeismicWavePC
# (seismicWavePCCases.jl) so that the six codes see exactly the same problem.
# Plain Julia (no packages): it is also loaded in the isolated SeismicWavePC environment.
#
# Cases
#   tilt0, tilt10, tilt20, tilt30, tilt45 : the flat point-force benchmark rotated by θ
#       (free surface n̂·r = 0, n̂ = (-sin θ, cos θ); source, force (along n̂) and
#       receivers rotate with it; components reported as (u_t, u_n) in the rotated frame)
#   full0 ... full45 : the same without free surface (grid-anisotropy control, flexOPT only)
#   hill : Gaussian hill z = H exp(-(x - xc)²/(2w²)) above z = 0, H = 3 km, w = 5 km,
#       xc = 0 (between the source at x = -10 km and the receivers); same source and
#       force as tilt0; receivers at x = -30 ... 10 km, "deep" at z = -10 km and
#       "shallow" 1 km below the local surface; components (u_t, u_n) = (u_x, u_z)

const GC = (vp=6000.0, vs=3464.0, rho=2700.0, halfWidth=60e3, voidCapCells=3, duration=9.0,
    sourceFlat=(-10e3, -20e3), receiversX=(-30e3, -20e3, -10e3, 0.0, 10e3),
    lines=(deep=-10e3, shallow=-1e3), force=1.0e10, sigma=500.0, cutoff=2.0,
    referenceDx=500.0, pointsPerSWavelength=16, rickerSharpness=1.5,
    hill=(height=3e3, width=5e3, center=0.0))

gcFrequency() = GC.rickerSharpness * GC.vs / (GC.pointsPerSWavelength * GC.referenceDx)
gcDelay() = 1.5 / gcFrequency()
"Source time function used by every code (SPECFEM2D/SeismicWavePC get it sampled)."
gcRicker(t) = (a = π * gcFrequency() * (t - gcDelay()); (1 - 2a^2) * exp(-a^2))

caseKind(case) = startswith(case, "tilt") ? :tilt : startswith(case, "full") ? :full :
    case == "hill" ? :hill : error("unknown case $case")
caseTheta(case) = caseKind(case) == :hill ? 0.0 :
    parse(Float64, replace(replace(case, "tilt" => ""), "full" => ""))
hasFreeSurface(case) = caseKind(case) != :full

"Unit tangent and normal of the rotated frame (θ = 0 for the hill)."
function caseFrame(case)
    θ = deg2rad(caseTheta(case))
    (t=(cos(θ), sin(θ)), n=(-sin(θ), cos(θ)))
end
"Flat-frame (x', z') → grid (x, z)."
function caseToGrid(case, xp, zp)
    f = caseFrame(case)
    (xp * f.t[1] + zp * f.n[1], xp * f.t[2] + zp * f.n[2])
end
"Free-surface height z_s(x) (material below); `nothing` without free surface."
function caseSurface(case)
    k = caseKind(case)
    k == :full && return nothing
    if k == :hill
        h = GC.hill
        return x -> h.height * exp(-(x - h.center)^2 / (2h.width^2))
    end
    θ = deg2rad(caseTheta(case))
    return x -> x * tan(θ)
end
"Level-set function φ(x, z), ≤ 0 in the material (exact for the plane: n̂·r)."
function casePhi(case)
    k = caseKind(case)
    k == :full && return (x, z) -> -1.0
    if k == :tilt
        n = caseFrame(case).n
        return (x, z) -> n[1] * x + n[2] * z
    end
    s = caseSurface(case)
    return (x, z) -> z - s(x)
end
caseSource(case) = caseToGrid(case, GC.sourceFlat...)
"Force direction (unit): along n̂ (vertical for the hill)."
caseForce(case) = caseFrame(case).n
"Receivers: Dict(:deep => [(x, z)...], :shallow => [...])."
function caseReceivers(case)
    if caseKind(case) == :hill
        s = caseSurface(case)
        return Dict(:deep => [(xr, GC.lines.deep) for xr in GC.receiversX],
            :shallow => [(xr, s(xr) + GC.lines.shallow) for xr in GC.receiversX])
    end
    Dict(line => [caseToGrid(case, xr, zl) for xr in GC.receiversX] for (line, zl) in pairs(GC.lines))
end
"Grid axes (multiples of dx) covering the case."
function caseGrid(case, dx)
    sx, sz = GC.sourceFlat
    W = GC.halfWidth
    if caseKind(case) == :hill
        x = collect(floor((sx - W) / dx) * dx:dx:ceil((sx + W) / dx) * dx)
        top = GC.hill.height + GC.voidCapCells * dx
        z = collect(floor((sz - W) / dx) * dx:dx:ceil(top / dx) * dx)
        return x, z
    end
    top = hasFreeSurface(case) ? GC.voidCapCells * dx : 20e3
    corners = [caseToGrid(case, a, b) for a in (sx - W, sx + W), b in (sz - W, top)]
    lo(d) = floor(minimum(c[d] for c in corners) / dx) * dx
    hi(d) = ceil(maximum(c[d] for c in corners) / dx) * dx
    collect(lo(1):dx:hi(1)), collect(lo(2):dx:hi(2))
end
"Gaussian source weights on the nodes of (x, z): [(i, k, weight)], Σ weight = 1."
function caseSourceWeights(case, x, z)
    s = caseSource(case)
    specs = [(i, k, exp(-((x[i] - s[1])^2 + (z[k] - s[2])^2) / (2GC.sigma^2)))
             for i in eachindex(x), k in eachindex(z) if hypot(x[i] - s[1], z[k] - s[2]) <= GC.cutoff * GC.sigma]
    total = sum(last, specs)
    [(i, k, w / total) for (i, k, w) in specs]
end
"(u_x, u_z) → (u_t, u_n) in the case frame."
function caseRotate(case, ux, uz)
    f = caseFrame(case)
    (ux .* f.t[1] .+ uz .* f.t[2], ux .* f.n[1] .+ uz .* f.n[2])
end

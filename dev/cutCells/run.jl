# dev/cutCells/run.jl -- checks and figures for src/cutCells/ (cut-cell
# volume integrals for free surfaces and internal discontinuities).
#
# Step `freeEnd1D` (stage 2): OPT3 rows of cut nodes on a 1-D bar whose free
# end falls between two nodes, against the exact solution (freeEnd1D.jl);
# cases eta0 eta025 eta05 eta075 (fraction of a cell beyond the last solid node).
#
# Step `geometry` (stage 1): exactness of the simplex rules, cut-cell moments
# against closed forms in 1-D/2-D/3-D, convergence of circle/sphere areas,
# tilted-plane geometry (normals, crossings, points+normals level set), and
# one figure per test interface. Every check asserts; the run ends with
# "ALL CHECKS PASSED".
#
# Loads only src/cutCells/ (in a small wrapper module), not the whole of
# flexOPT, so it takes seconds. src/flexOPT.jl includes the same files.
#
# Run:
#   julia --project=. dev/cutCells/run.jl geometry
#   julia --project=. dev/cutCells/run.jl freeEnd1D            # all η cases
#   julia --project=. dev/cutCells/run.jl freeEnd1D eta05
# REPL:
#   julia> include("dev/cutCells/run.jl"); cutCells("geometry")

if !@isdefined(flexopt_root)
    import Pkg
    function find_flexopt_root(start_dir=@__DIR__)
        directory = abspath(start_dir)
        while true
            isfile(joinpath(directory, "src", "flexOPT.jl")) && return directory
            parent = dirname(directory)
            parent == directory && error("Cannot locate flexOPT; set ENV[\"FLEXOPT_ROOT\"]")
            directory = parent
        end
    end
    global flexopt_root = haskey(ENV, "FLEXOPT_ROOT") ? abspath(expanduser(ENV["FLEXOPT_ROOT"])) : find_flexopt_root()
    Pkg.activate(flexopt_root)
end
include(joinpath(flexopt_root, "dev", "common", "devBootstrap.jl"))

cutCells(args...) = (global APP_ARGS = collect(String, string.(args)); include(joinpath(@__DIR__, "run.jl")))

DEV = parseDevArgs(@isdefined(APP_ARGS) ? APP_ARGS : ARGS; steps=["geometry", "freeEnd1D"], defaultStep="geometry",
    cases=["eta0", "eta025", "eta05", "eta075"])
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")

if !@isdefined(CutCellsCheck) || DEV.reload
    @eval module CutCellsCheck
        const CUTCELLS_GEOMETRY_ONLY = true
        include(joinpath(Main.flexopt_root, "src", "cutCells", "cutCells.jl"))
    end
end
using .CutCellsCheck
using StaticArrays, LinearAlgebra, SparseArrays, Printf, CairoMakie
CairoMakie.activate!(type="png")

# --------------------------------------------------------------- utilities
"Gauss–Legendre nodes/weights on [a, b] (Golub–Welsch)."
function gaussLegendre(n, a, b)
    β = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    x, w = E.values, 2 .* E.vectors[1, :] .^ 2
    ((b - a) / 2 .* x .+ (a + b) / 2, (b - a) / 2 .* w)
end
"∫_lo^hi f with the interval split at `breaks` (12-point rule per piece)."
function piecewiseQuad(f, lo, hi, breaks)
    cuts = sort(unique(vcat(lo, hi, filter(b -> lo < b < hi, breaks))))
    total = 0.0
    for k in 1:length(cuts)-1
        x, w = gaussLegendre(12, cuts[k], cuts[k+1])
        total += sum(w .* f.(x))
    end
    total
end
function check(name, value, reference; atol=1e-12)
    err = maximum(abs, value .- reference)
    @printf("  %-58s max error %.2e  (tol %.0e)\n", name, err, atol)
    err <= atol || error("check failed: $name")
end
cornerArray(f, N) = [f(Tuple(k)) for k in CartesianIndices(ntuple(_ -> 0:1, N))]

# ================================================================ geometry
if "geometry" in DEV.steps
    DEV_STEP[] = "geometry"
    println("\n== simplex rules")
    tri = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(0.0, 1.0)]
    m = zeros(6, 6)
    addSimplexMonomialMoments!(m, tri, 5)
    check("triangle x^a z^b = a!b!/(a+b+2)!, a,b<=5", m,
        [factorial(a) * factorial(b) / factorial(a + b + 2) for a in 0:5, b in 0:5])
    tet = [SVector(0.0, 0, 0), SVector(1.0, 0, 0), SVector(0.0, 1, 0), SVector(0.0, 0, 1)]
    m = zeros(4, 4, 4)
    addSimplexMonomialMoments!(m, tet, 3)
    check("tetrahedron x^a y^b z^c, a,b,c<=3", m,
        [factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + 3) for a in 0:3, b in 0:3, c in 0:3])

    println("\n== cell moments")
    for N in 1:3
        check("uncut $(N)-D cell, Π 1/(a+1)", cellMoments(cornerArray(_ -> -1.0, N), 3), fullCellMoments(N, 3))
        check("empty $(N)-D cell", cellMoments(cornerArray(_ -> 1.0, N), 3), zeros(ntuple(_ -> 4, N)))
    end
    for ξ in (0.0, 0.3, 0.71, 1.0)
        check("1-D, material y <= $ξ", cellMoments(cornerArray(k -> k[1] - ξ, 1), 5),
            [ξ^(a + 1) / (a + 1) for a in 0:5])
    end
    # 2-D: material z <= h(x) = z0 + s x (φ = z - h is linear, so exact)
    for (z0, s) in ((0.3, 0.4), (0.2, 1.5), (0.8, -0.6), (-0.2, 0.9), (0.5, 0.0))
        h(x) = z0 + s * x
        breaks = s == 0 ? Float64[] : [(0 - z0) / s, (1 - z0) / s]
        ref = [piecewiseQuad(x -> x^a * clamp(h(x), 0, 1)^(b + 1) / (b + 1), 0.0, 1.0, breaks) for a in 0:3, b in 0:3]
        mm = cellMoments(cornerArray(k -> k[2] - h(k[1]), 2), 3)
        check("2-D, material z <= $z0 + $s x", mm, ref)
        check("2-D, '−' + '+' = whole cell ($z0, $s)",
            mm .+ cellMoments(cornerArray(k -> k[2] - h(k[1]), 2), 3; side=:plus), fullCellMoments(2, 3))
    end
    # 3-D: material z <= c0 + c1 x + c2 y, plane inside the cell (prism cases occur)
    let c0 = 0.3, c1 = 0.2, c2 = 0.3
        h(x, y) = c0 + c1 * x + c2 * y
        xg, wg = gaussLegendre(8, 0.0, 1.0)
        ref = [sum(wg[i] * wg[j] * xg[i]^a * xg[j]^b * h(xg[i], xg[j])^(c + 1) / (c + 1)
            for i in 1:8, j in 1:8) for a in 0:2, b in 0:2, c in 0:2]
        check("3-D, material z <= 0.3 + 0.2x + 0.3y", cellMoments(cornerArray(k -> k[3] - h(k[1], k[2]), 3), 2), ref)
    end
    let n = normalize(SVector(0.3, -0.5, 0.8)), c = 0.1
        φc = cornerArray(k -> dot(n, SVector{3,Float64}(k)) - c, 3)
        check("3-D oblique plane, '−' + '+' = whole cell",
            cellMoments(φc, 3) .+ cellMoments(φc, 3; side=:plus), fullCellMoments(3, 3))
    end

    println("\n== convergence: circle (2-D) and sphere (3-D), φ = |r| - R")
    function enclosed(N, R, h)
        axis = collect(-0.5:h:0.5)
        ls = levelSetFromFunction(r -> norm(r) - R, ntuple(_ -> axis, N))
        geom = CutCellGeometry(ls; maxDegree=1)
        sum(geom.fraction) * h^N
    end
    for (N, R, hs, exact) in ((2, 0.37, (1/8, 1/16, 1/32, 1/64, 1/128), π * 0.37^2),
                              (3, 0.37, (1/8, 1/16, 1/32), 4π / 3 * 0.37^3))
        errors = [abs(enclosed(N, R, h) - exact) for h in hs]
        rates = [log2(errors[k] / errors[k+1]) for k in 1:length(hs)-1]
        for (h, e) in zip(hs, errors)
            @printf("  %d-D  h = 1/%-4d  |error| = %.3e\n", N, round(Int, 1 / h), e)
        end
        @printf("  %d-D  observed orders %s\n", N, join([@sprintf("%.2f", r) for r in rates], ", "))
        rates[end] > 1.7 || error("check failed: $(N)-D enclosed volume does not converge at second order")
    end

    println("\n== tilted plane (2-D), 20°")
    θ = deg2rad(20.0)
    n̂ = SVector(-sin(θ), cos(θ))
    hgrid = 0.5
    x = collect(-10.0:hgrid:10.0)
    z = collect(-8.0:hgrid:6.0)
    ls = levelSetFromFunction(r -> dot(n̂, r), (x, z))
    geom = CutCellGeometry(ls; maxDegree=4)
    nodes = cutNodes(geom)
    solid = sideMask(ls)
    @printf("  %d cut cells, %d cut nodes (%d solid, %d void)\n", length(geom.moments), length(nodes),
        count(I -> solid[I], nodes), count(I -> !solid[I], nodes))
    crossings = reduce(vcat, [stencilCrossings(ls, I) for I in nodes])
    check("normals at all crossings = (-sin θ, cos θ)", maximum(c -> norm(c.normal - n̂), crossings), 0.0; atol=1e-12)
    check("crossing points on the plane", maximum(c -> abs(dot(n̂, c.point)), crossings), 0.0; atol=1e-12)
    # area of {n̂·r <= 0} ∩ box: polygon clipping of the box
    function clippedBoxArea(n, x0, x1, z0, z1)
        poly = [SVector(x0, z0), SVector(x1, z0), SVector(x1, z1), SVector(x0, z1)]
        out = SVector{2,Float64}[]
        for k in eachindex(poly)
            p, q = poly[k], poly[mod1(k + 1, length(poly))]
            fp, fq = dot(n, p), dot(n, q)
            fp <= 0 && push!(out, p)
            (fp <= 0) != (fq <= 0) && push!(out, p + fp / (fp - fq) * (q - p))
        end
        0.5 * abs(sum(out[k][1] * out[mod1(k + 1, end)][2] - out[mod1(k + 1, end)][1] * out[k][2] for k in eachindex(out)))
    end
    check("Σ fractions · h² = area below the plane", sum(geom.fraction) * hgrid^2,
        clippedBoxArea(n̂, first(x), last(x), first(z), last(z)); atol=1e-9)
    samples = [SVector(t * cos(θ), t * sin(θ)) for t in -15.0:0.25:15.0]
    lsPN = levelSetFromPointsNormals(samples, fill(n̂, length(samples)), (x, z))
    check("points + normals level set = n̂·r for a plane", lsPN.φ, ls.φ; atol=1e-12)
    probe = SVector(1.3, 0.2)
    check("multilinear φ at an off-node point (plane: exact)", levelSetValue(ls, probe), dot(n̂, probe); atol=1e-12)

    println("\n== small cut cells: void nodes merged into a solid neighbour")
    aggregation = sliverAggregation(geom)
    voidCut = [I for I in nodes if !solid[I]]
    @printf("  %d void cut nodes, %d merged\n", length(voidCut), length(aggregation))
    length(aggregation) == length(voidCut) || error("check failed: every void cut node must be merged")
    for agg in aggregation
        solid[agg.target] || error("check failed: target $(agg.target) is not solid")
        all(solid[node] for (node, _) in agg.extension) || error("check failed: extension uses a void node")
        o = SVector{2,Float64}(Tuple(agg.target - agg.node) .* hgrid)
        dot(normalize(o), -n̂) > 0 || error("check failed: merge direction points away from the solid")
    end
    println("  targets solid, extensions solid, directions inward: ok")
    let A = sparse([2.0 -1 0; -1 2 -1; 0 -1 2]), L = sparse([1.0 0 0; 0 1 0; 0 0 1]), R = copy(L)
        prepared = (A_unknown=A, L_known=L, R_force=R, NField=1, NpointsSpace=3, spaceShape=(3,))
        out = applyAggregation(prepared, [(node=CartesianIndex(3), target=CartesianIndex(2), extension=[(CartesianIndex(2), 2.0), (CartesianIndex(1), -1.0)])])
        check("applyAggregation on a 3-node system (A)", Matrix(out.A_unknown), [2.0 -1 0; -1 1 1; 1 -2 1])
        check("applyAggregation on a 3-node system (L)", Matrix(out.L_known), [1.0 0 0; 0 1 1; 0 0 0])
    end

    println("\n== figures")
    function geometryFigure(ls, geom, title; arrowsEvery=3)
        xa, za = ls.axes
        nodes = cutNodes(geom)
        solid = sideMask(ls)
        fig = Figure(size=(1100, 520))
        ax = Axis(fig[1, 1]; title=title, xlabel="x", ylabel="z", aspect=DataAspect())
        xc = (xa[1:end-1] .+ xa[2:end]) ./ 2
        zc = (za[1:end-1] .+ za[2:end]) ./ 2
        hm = heatmap!(ax, xc, zc, geom.fraction; colormap=:grays, colorrange=(0, 1))
        Colorbar(fig[1, 2], hm; label="solid fraction of the cell")
        contour!(ax, xa, za, ls.φ; levels=[0.0], color=:dodgerblue, linewidth=2)
        sn = [I for I in nodes if solid[I]]
        vn = [I for I in nodes if !solid[I]]
        scatter!(ax, [xa[I[1]] for I in sn], [za[I[2]] for I in sn]; color=:crimson, markersize=7, label="cut node, φ ≤ 0")
        scatter!(ax, [xa[I[1]] for I in vn], [za[I[2]] for I in vn]; color=:orange, marker=:utriangle, markersize=7, label="cut node, φ > 0")
        cr = reduce(vcat, [stencilCrossings(ls, I) for I in sn[1:arrowsEvery:end]]; init=NamedTuple[])
        if !isempty(cr)
            arrows!(ax, [c.point[1] for c in cr], [c.point[2] for c in cr],
                [c.normal[1] for c in cr], [c.normal[2] for c in cr]; lengthscale=0.6, color=:dodgerblue)
        end
        axislegend(ax; position=:lb)
        fig
    end
    devSave(geometryFigure(ls, geom, "tilted plane 20°: cut cells, cut nodes, normals at crossings"), "tilted_plane_20deg")
    xh = collect(-10.0:0.5:10.0)
    zh = collect(-6.0:0.5:6.0)
    lsHill = levelSetFromHeight(x -> 3.0 * exp(-x^2 / (2 * 2.5^2)), (xh, zh))
    devSave(geometryFigure(lsHill, CutCellGeometry(lsHill), "Gaussian hill (σ = 2.5, height 3): cut cells, cut nodes, normals"), "gaussian_hill")

    println("\nALL CHECKS PASSED")
end

# =============================================================== freeEnd1D
if "freeEnd1D" in DEV.steps
    DEV_STEP[] = "freeEnd1D"
    include(joinpath(@__DIR__, "freeEnd1D.jl"))
end

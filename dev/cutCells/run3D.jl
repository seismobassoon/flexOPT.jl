# dev/cutCells/run3D.jl -- flexOPT 3-D topography runs with cut cells (or zigzag)
#
# Cases (geometryCases3D.jl): flat3D, oblique10, oblique20 (planes; Gaussian hill and
# Kirishima come next). Same scheme switches as the 2-D runs (env.jl):
#   FLEXOPT_CUT_SCHEME=OPT3|OPT5, FLEXOPT_CUT_GEOMETRY=cut|zigzag, FLEXOPT_CPU_RECIPE=1.
#
# Numerics
#  * OPT3 recipe of famousEquations("3DsismoTimeIsoHeteroForce"), weak form, Δ = (1,1,1,1).
#  * Operators on the model box + Cerjan padding (4 sides and the bottom, none on top), as
#    dev/greenland_nov2026. The level set, cut cells and cut rows are built on the WHOLE
#    (padded) grid; the Cerjan weight exp(-d·dist²) on past-time columns of padding points
#    is re-applied to the installed cut rows (cutCellRows does not know about Cerjan).
#  * Void nodes: extension rows (u = Lagrange extrapolation), sliver merging, far void
#    pinned -- exactly as in 2-D.
#  * Time marching: predictor-corrector (no LU in 3-D), with the lumped-future blocks of
#    flexOPT except at constraint points (extended and pinned void nodes), where the block
#    is the diagonal coefficient only: the lumped block of an extrapolation row is
#    1 - Σw = 0. Only the receiver nodes are stored.
#
# Run (flexOPT root):
#   FLEXOPT_CPU_RECIPE=1 julia --project=. -t auto dev/cutCells/run3D.jl flat3D oblique10 oblique20
#   options: --redo (recompute), --small (5 km box, 6 Cerjan cells, 2 s), --tiny (2 km, 0.6 s): pipeline tests
# Cache: data/dev/cutCells/run3D_<hash>.jld2; traces "<line>_<i>" => (u1, u2, un).

if !@isdefined(flexopt_root)
    import Pkg
    global flexopt_root = abspath(joinpath(@__DIR__, "..", ".."))
    Pkg.activate(flexopt_root)
end
include(joinpath(@__DIR__, "env.jl"))
include(joinpath(@__DIR__, "geometryCases3D.jl"))
using StaticArrays, SparseArrays, LinearAlgebra, Printf, JLD2, DrWatson
using .flexOPT: CerjanBoundarySpec, BoundaryConditionSet

const RUN3D_VERSION = "v1-2026-11-26"
SMALL3D = "--small" in ARGS
TINY3D = "--tiny" in ARGS
REDO3D = "--redo" in ARGS
cases3D = filter(a -> !startswith(a, "--"), ARGS)
isempty(cases3D) && (cases3D = collect(CASES3D))
size3D = TINY3D ? (halfWidth=2e3, below=2.5e3, padding=3, duration=0.6) :
    SMALL3D ? (halfWidth=5e3, below=5e3, padding=6, duration=2.0) :
    (halfWidth=GC3.halfWidth, below=GC3.below, padding=GC3.cerjanCells, duration=GC3.duration)
haskey(ENV, "FLEXOPT_3D_DURATION") && (size3D = merge(size3D, (duration=parse(Float64, ENV["FLEXOPT_3D_DURATION"]),)))
pc3D = (tolerance=1e-10, maxIterations=60, relaxation=1.0)

# ------------------------------------------------------------------ recipe and cut-cell scheme
optParameters3D = Dict{String,Any}(
    "famousEquationType" => "3DsismoTimeIsoHeteroForce",
    "Δ" => (1.0, 1.0, 1.0, 1.0), "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3, "supplementaryOrder" => 2,
    "variationalForm" => DEFAULT_FORM[],
    "taylor_inverse_mode" => :weak_operator_optimized,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "recipe_backend" => backend,
)
optParameters3D = schemeParameters(optParameters3D, CUT_SCHEME)
recipe3D = cachedRecipe(optParameters3D, "elastic3D_$(CUT_SCHEME)")
@time scheme3D = flexOPT.CutCellScheme(optParameters3D)

"Lumped-future inverse blocks; at `diagonalPoints` only the point's own coefficients."
function predictorBlocks(A::SparseMatrixCSC, NF::Int, NP::Int, diagonalPoints::BitVector)
    blocks = zeros(NF, NF, NP)
    rv, nzv = rowvals(A), nonzeros(A)
    for c in 1:size(A, 2)
        unknownField = (c - 1) ÷ NP + 1
        colPoint = (c - 1) % NP + 1
        for k in nzrange(A, c)
            r = rv[k]
            point = (r - 1) ÷ NF + 1
            (diagonalPoints[point] && point != colPoint) && continue
            blocks[(r - 1) % NF + 1, unknownField, point] += nzv[k]
        end
    end
    inverse = similar(blocks)
    worst = 0.0
    for p in 1:NP
        B = SMatrix{3,3}(@view blocks[:, :, p])
        inverse[:, :, p] .= inv(B)
        worst = max(worst, cond(B))
    end
    (npoints=NP, nfield=NF, inverse_blocks=inverse, worstCondition=worst)
end

"""
    bicgstabLumped!(x, A, b, P, w; tol, maxit)

Right-preconditioned BiCGStab for A x = b, preconditioner = the lumped-future blocks P
(rows point-major, unknowns field-major, as flexOPT). `x` holds the initial guess (the
previous step). Returns (iterations, final relative residual).
"""
function bicgstabLumped!(x, A, b, P, w; tol=1e-10, maxit=60)
    nb = max(norm(b), eps())
    mul!(w.r, A, x); w.r .= b .- w.r
    rel = norm(w.r) / nb
    rel <= tol && return 0, rel
    w.rhat .= w.r; w.p .= 0; w.v .= 0
    ρ = α = ω = 1.0
    for it in 1:maxit
        ρnew = dot(w.rhat, w.r)
        β = (ρnew / ρ) * (α / ω); ρ = ρnew
        w.p .= w.r .+ β .* (w.p .- ω .* w.v)
        flexOPT.applyLumpedFuturePredictor!(w.y, P, w.p)
        mul!(w.v, A, w.y)
        α = ρ / dot(w.rhat, w.v)
        w.s .= w.r .- α .* w.v
        if norm(w.s) / nb <= tol
            x .+= α .* w.y
            return it, norm(w.s) / nb
        end
        flexOPT.applyLumpedFuturePredictor!(w.z, P, w.s)
        mul!(w.t, A, w.z)
        ω = dot(w.t, w.s) / dot(w.t, w.t)
        x .+= α .* w.y .+ ω .* w.z
        w.r .= w.s .- ω .* w.t
        rel = norm(w.r) / nb
        rel <= tol && return it, rel
    end
    return maxit, rel
end

"Trilinear weights and the 8 whole-grid nodes around p."
function trilinear(xw, yw, zw, p)
    h = xw[2] - xw[1]
    f = ((p[1] - xw[1]) / h, (p[2] - yw[1]) / h, (p[3] - zw[1]) / h)
    i0 = floor.(Int, f) .+ 1
    a = f .- (i0 .- 1)
    nodes = CartesianIndex{3}[]; w = Float64[]
    for di in 0:1, dj in 0:1, dk in 0:1
        push!(nodes, CartesianIndex(i0[1] + di, i0[2] + dj, i0[3] + dk))
        push!(w, (di == 1 ? a[1] : 1 - a[1]) * (dj == 1 ? a[2] : 1 - a[2]) * (dk == 1 ? a[3] : 1 - a[3]))
    end
    nodes, w
end

function run3D(case)
    dx = GC3.dx
    dt = 0.20 * dx / (sqrt(3) * GC3.vp * 2)
    g = case3DGrid(case; halfWidth=size3D.halfWidth, below=size3D.below, padding=size3D.padding)
    Nmodel = (length(g.x), length(g.y), length(g.z))
    Nwhole = (length(g.xw), length(g.yw), length(g.zw))
    P = size3D.padding
    padLower = CartesianIndex(P, P, P)
    @printf("  model %s, whole %s (%.2f M nodes, %.2f M unknowns), Δx = %g m, Δt = %.5f s\n",
        string(Nmodel), string(Nwhole), prod(Nwhole) / 1e6, 3prod(Nwhole) / 1e6, dx, dt)

    # material (homogeneous, scaled), operators with Cerjan padding
    scale2 = (dt / dx)^2
    mk(v, n) = fill(v, n)
    models = [mk(1.0, Nmodel), mk((GC3.vp^2 - 2GC3.vs^2) * scale2, Nmodel), mk(GC3.vs^2 * scale2, Nmodel)]
    modelsWhole = [mk(1.0, Nwhole), mk((GC3.vp^2 - 2GC3.vs^2) * scale2, Nwhole), mk(GC3.vs^2 * scale2, Nwhole)]
    damping = 3.0 / P^2
    cerjan = CerjanBoundarySpec((P, P, P), (P, P, 0); damping=damping)
    timeMarching = recipe3D["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching
    family = (models=models, modelPoints=getModelPoints(models[1], 3, timeMarching),
        Δ=(1.0, 1.0, 1.0, 1.0), modelName="cutcells3d_$(case)_$(TINY3D ? "tiny" : SMALL3D ? "small" : "std")")
    tBuild = @elapsed numerical = numericalOperatorConstruction(Dict{String,Any}(
        "optRec" => recipe3D, "modelFam" => family, "absorbingBoundaries" => nothing,
        "maskedRegionInSpace" => nothing,
        "boundaryConditions" => BoundaryConditionSet(free_surface=nothing, cerjan=cerjan, material_mask=nothing),
        "representation" => "matrixfree"))["numOperators"]
    volume = prepareLinearSystem(numerical)
    numerical = nothing
    volume.spaceShape == Nwhole || error("operator grid $(volume.spaceShape) differs from $Nwhole")
    NP, NF = volume.NpointsSpace, volume.NField
    @printf("  volume operator: %.0f s, nnz A %.1f M, L %.1f M, R %.1f M\n", tBuild,
        nnz(volume.A_unknown) / 1e6, nnz(volume.L_known) / 1e6, nnz(volume.R_force) / 1e6)

    # geometry on the whole grid
    φ = case3DPhi(case)
    ls = flexOPT.levelSetFromFunction(r -> φ(r[1], r[2], r[3]), (g.xw, g.yw, g.zw))
    tGeom = @elapsed geom = flexOPT.CutCellGeometry(ls; maxDegree=scheme3D.maxDegree, zigzag=(CUT_GEOMETRY == "zigzag"))
    solid = flexOPT.sideMask(geom.levelSet)
    cut = flexOPT.cutNodes(geom)
    cutSet = Set(cut)
    aggregation = flexOPT.sliverAggregation(geom; extension=CUT_EXTENSION, reach=CUT_REACH, maxSteps=CUT_REACH + 3)
    extended = Set(a.node for a in aggregation)
    farVoid = [I for I in CartesianIndices(Nwhole) if !solid[I] && !(I in cutSet) && !(I in extended)]
    nMerged = count(a -> a.aggregate, aggregation)
    @printf("  geometry %.0f s: %d cut nodes (%d solid), %d void merged, %d extended, %d pinned\n", tGeom,
        length(cut), count(I -> solid[I], cut), nMerged, length(aggregation) - nMerged, length(farVoid))

    # self-check at the source node, cut rows, Cerjan weights on their past-time columns
    srcWhole = CartesianIndex(argmin(abs.(g.xw .- case3DSource(case)[1])), argmin(abs.(g.yw .- case3DSource(case)[2])),
        argmin(abs.(g.zw .- case3DSource(case)[3])))
    selfCheck = flexOPT.cutCellSelfCheck(scheme3D, volume, modelsWhole, srcWhole)
    @show selfCheck
    maximum(values(selfCheck)) < selfCheckTolerance(recipe3D) || error("cut-cell rows of an uncut node differ from the volume rows: $selfCheck")
    println("  computing cut-cell rows for $(length(cut)) nodes ..."); flush(stdout)
    tRows = @elapsed rows = flexOPT.cutCellRows(scheme3D, geom, cut, modelsWhole, volume)
    CIw = CartesianIndices(Nwhole)
    one3, Nm = CartesianIndex(1, 1, 1), CartesianIndex(Nmodel)
    nWeighted = 0
    for (row, entries) in rows.L
        for k in eachindex(entries)
            col, v = entries[k]
            Im = CIw[(col - 1) % NP + 1] - padLower
            d2 = distance2_point_to_box(Im, one3, Nm)
            d2 > 0 && (entries[k] = (col, v * exp(-damping * d2)); nWeighted += 1)
        end
    end
    @printf("  cut-cell rows %.0f s (%d L entries Cerjan-weighted)\n", tRows, nWeighted)
    @printf("  memory before installing rows: %.1f GB (max RSS)\n", Sys.maxrss() / 1e9)
    system = flexOPT.pinRows(flexOPT.applyAggregation(flexOPT.installRows(volume, rows), aggregation), farVoid)
    volume = nothing; rows = nothing
    GC.gc()

    # condensation of the void constraints: u = E u_r, u_r = unknowns of the material nodes;
    # extended void nodes u_v = Σ w u_J (J material), far void u = 0
    LIw = LinearIndices(Nwhole)
    solidPts = [LIw[I] for I in CartesianIndices(Nwhole) if solid[I]]
    Ns = length(solidPts)
    reduced = zeros(Int, NP); reduced[solidPts] .= 1:Ns
    Ei, Ej, Ev = Int[], Int[], Float64[]
    for f in 1:NF, (k, p) in enumerate(solidPts)
        push!(Ei, p + (f - 1) * NP); push!(Ej, k + (f - 1) * Ns); push!(Ev, 1.0)
    end
    for a in aggregation, f in 1:NF, (J, w) in a.extension
        kJ = reduced[LIw[J]]
        kJ > 0 || error("extension of $(a.node) reads a void node $J")
        push!(Ei, LIw[a.node] + (f - 1) * NP); push!(Ej, kJ + (f - 1) * Ns); push!(Ev, w)
    end
    E = sparse(Ei, Ej, Ev, NP * NF, Ns * NF)
    rowsS = [(p - 1) * NF + e for p in solidPts for e in 1:NF]
    tCond = @elapsed begin
        Ar = system.A_unknown[rowsS, :] * E
        Lr = system.L_known[rowsS, :]
    end
    tPred = @elapsed predictor = predictorBlocks(Ar, NF, Ns, falses(Ns))
    @printf("  condensed system: %d material nodes (%.0f s); predictor %.0f s, worst block condition %.2e\n",
        Ns, tCond, tPred, predictor.worstCondition)

    # source: Gaussian weights (model-box nodes) → whole grid, 9 force columns per node
    specs = case3DSourceWeights(case, g.x, g.y, g.z)
    n = case3DForce(case)
    NFP, NFF = system.NforcePoints, system.NForceField
    srcCols = Int[]; srcVals = Float64[]           # value = weight · direction component; × wavelet(t)
    for (i, j, k, w) in specs
        Iw = CartesianIndex(i, j, k) + padLower
        solid[Iw] || error("source touches the void")
        lp = LIw[Iw]
        for s in 1:3, c in 1:3
            push!(srcCols, lp + (c - 1) * NFP + (s - 1) * NFP * NFF)
            push!(srcVals, w * n[c])
        end
    end
    Rsrc = system.R_force[rowsS, srcCols]
    system = nothing; GC.gc()
    sourceScale = GC3.force * dt^2 / (GC3.rho * dx^3)
    slotOf = repeat(repeat(1:3, inner=3), length(specs))     # s for each column (c fastest)

    # receivers
    receivers = case3DReceivers(case)
    recNodes = Dict{String,Tuple{Vector{Int},Vector{Float64}}}()
    for line in (:deep, :shallow), (ir, p) in enumerate(receivers[line])
        nodes, w = trilinear(g.xw, g.yw, g.zw, p.r)
        if !all(I -> checkbounds(Bool, solid, I), nodes)
            (TINY3D || SMALL3D) || error("receiver $line $ir outside the grid")
            continue                                   # test sizes: receivers off the small box are skipped
        end
        all(solid[I] for I in nodes) || error("receiver $line $ir interpolates from the void")
        recNodes["$(line)_$(ir)"] = ([LIw[I] for I in nodes], w)
    end

    # march: A_r u_r = -L_r u_known + R_r f, BiCGStab preconditioned by the lumped blocks
    steps = ceil(Int, size3D.duration / dt)
    known = zeros(NP, NF, 2)
    nr = Ns * NF
    ur = zeros(nr); rhs = zeros(nr); u = zeros(NP * NF)
    work = (r=zeros(nr), rhat=zeros(nr), p=zeros(nr), v=zeros(nr), s=zeros(nr), t=zeros(nr), y=zeros(nr), z=zeros(nr))
    fvals = zeros(length(srcCols))
    rec = Dict(k => zeros(3, steps + 1) for k in keys(recNodes))
    iters = zeros(Int, steps); resid = zeros(steps)
    t0 = time(); stopped = false; lastStep = steps
    for nstep in 1:steps
        for q in eachindex(fvals)
            fvals[q] = sourceScale * srcVals[q] * gc3Ricker((nstep - 3 + slotOf[q]) * dt)
        end
        mul!(rhs, Lr, vec(known), -1.0, 0.0)
        mul!(rhs, Rsrc, fvals, 1.0, 1.0)
        used, rel = bicgstabLumped!(ur, Ar, rhs, predictor, work; tol=pc3D.tolerance, maxit=pc3D.maxIterations)
        iters[nstep] = used; resid[nstep] = rel
        mul!(u, E, ur)
        amp = maximum(abs, u)
        if !isfinite(amp) || amp > 1e8
            @warn "stopping: blow-up" nstep amp
            stopped = true; lastStep = nstep - 1; break
        end
        known[:, :, 1] .= @view known[:, :, 2]
        known[:, :, 2] .= reshape(u, NP, NF)
        for (key, (lps, w)) in recNodes, c in 1:3
            rec[key][c, nstep+1] = sum(w[q] * u[lps[q] + (c - 1) * NP] for q in eachindex(lps))
        end
        if nstep % max(1, steps ÷ 20) == 0
            el = time() - t0
            @printf("  step %5d / %d  t = %.2f s  max|u| %.2e  BiCGStab it %d res %.1e  (%.0f s, ETA %.0f s)\n",
                nstep, steps, nstep * dt, amp, used, resid[nstep], el, el / nstep * (steps - nstep))
        end
    end
    stopped && error("propagation stopped early (blow-up) for $case")
    marchTime = time() - t0
    traces = Dict{String,Any}()
    for line in (:deep, :shallow), (ir, p) in enumerate(receivers[line])
        haskey(rec, "$(line)_$(ir)") || continue
        r = rec["$(line)_$(ir)"]
        u1, u2, un = case3DRotate(case, r[1, :], r[2, :], r[3, :])
        traces["$(line)_$(ir)"] = (xp_km=p.xp / 1e3, yp_km=p.yp / 1e3, line=String(line), u1=u1, u2=u2, un=un)
    end
    Dict{String,Any}("case" => case, "times" => collect(0:steps) .* dt, "traces" => traces,
        "selfCheck" => selfCheck, "nCut" => length(cut), "nMerged" => nMerged,
        "nExtended" => length(aggregation) - nMerged, "solverIterations" => iters, "solverResiduals" => resid,
        "marchTime" => marchTime, "rowsTime" => tRows, "wholeShape" => Nwhole)
end

params3D(case) = Dict{String,Any}("stage" => "run3D", "case" => case, "scheme" => CUT_SCHEME,
    "geometry" => CUT_GEOMETRY, "extension" => String(CUT_EXTENSION), "dx" => GC3.dx,
    "halfWidth" => size3D.halfWidth, "below" => size3D.below, "padding" => size3D.padding,
    "duration" => size3D.duration, "codeVersion" => RUN3D_VERSION)
cacheFile3D(case) = datadir("dev/cutCells", "run3D_$(commonBatchs.hash_parameters(params3D(case))).jld2")

results3D = Dict{String,Any}()
for case in cases3D
    println("\n== 3-D $(CUT_VARIANT): $case (θ = $(case3DTheta(case))°, φ = $(case3DAzimuth(case))°)",
        TINY3D ? "  [--tiny]" : SMALL3D ? "  [--small]" : "")
    REDO3D && isfile(cacheFile3D(case)) && rm(cacheFile3D(case))
    results3D[case] = myProduceOrLoad(cfg -> run3D(case), params3D(case), "dev/cutCells", "run3D")
    tr = results3D[case]["traces"]["deep_3"]
    @printf("  check above the source: max(|u_1|,|u_2|)/max|u_n| = %.3f (≈ 0)\n",
        max(maximum(abs, tr.u1), maximum(abs, tr.u2)) / max(maximum(abs, tr.un), eps()))
    println("  -> ", cacheFile3D(case))
end

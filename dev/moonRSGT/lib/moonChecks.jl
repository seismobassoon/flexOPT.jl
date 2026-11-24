# dev/moonRSGT/lib/moonChecks.jl -- validation of one catalogue model
#
# For a few sources (on the depth grid, off-grid, shallow) and a fixed
# station, compare the catalogue synthetics (receiver-side disk, Kᵀ,
# interpolated in depth AND distance) with an independent source-side
# computation (`computeSGTDisk(...; family = :moment)` at the exact source
# radius, exact distance, no interpolation, solved with K).  The difference
# measures interpolation + reciprocity together.

using Printf

"""
    moonCheck(ndFile, catFile, params; recipes, configs, outDir, station, M, sources)

Prints the relative L2 differences (spectra and Gaussian-filtered
displacement, per component) and saves a wiggle overlay figure
(needs CairoMakie loaded in Main).  Returns the table of errors.
"""
function moonCheck(ndFile, catFile, params; recipes, configs, outDir,
                   station=(-3.0, -23.4),
                   M=doubleCouple(strike=40, dip=55, rake=-70, M0=1e15),
                   sources=[(name="on grid, 300 km", src=(20.0, 15.0, 300.0)),
                            (name="off grid, 913 km", src=(-25.0, 70.0, 913.0)),
                            (name="off grid, 9 km (crust)", src=(10.0, -60.0, 9.0))],
                   comps=(:Z, :R, :T), halfDuration=10.0)
    cat = loadMoonRSGT(catFile)
    mm = moonModelFromND(ndFile; fitTol=params["fitTol"])
    psv, sh = moonSystems(mm.model, params; recipes=recipes, configs=configs)
    R = mm.info.R
    Mcode = collect(Float64, M) ./ 1e18
    Rcv = pointFromLatLon(R, 0.0, station...)
    rows = NamedTuple[]
    traces = []
    for s in sources
        lat, lon, dep = s.src
        Q = pointFromLatLon(R, dep, lat, lon)
        Rs, _ = sourceToPole(:moment, Mcode, Q)
        θp, _, _ = frameOf(Rs, Rcv)
        t0 = time()
        dS = computeSGTDisk(psv, sh, mm.model; family=:moment, rs_km=R - dep,
                            radii=(r_km=[R], side=[:upper]), θdeg=[θp],
                            tlen=params["tlen"], fmax=params["fmax"], omegaiTlen=params["omegaiTlen"],
                            p=sourceOrderOf(params["scheme"]), lmaxOf=moonLmaxOf(params, R), lcap=5000,
                            verbose=false)
        tref = time() - t0
        cs = moonSpectra(cat; station=station, source=s.src, M=M, comps=comps)
        for c in comps
            n = collect(receiverDirection(c, Q, Rcv))
            ref = [1e3 * dot(n, sgtAt(dS, k, Q, Mcode, Rcv; which=:disp)) for k in eachindex(dS.ω)]
            eS = norm(cs.spec[c] .- ref) / max(norm(ref), 1e-300)
            t, uc = moonTimeSeries(cs.ω, cs.ωI, cs.tlen, cs.spec[c]; halfDuration=halfDuration)
            _, ur = moonTimeSeries(dS.ω, dS.ωI, dS.tlen, ref; halfDuration=halfDuration)
            eT = norm(uc .- ur) / max(norm(ur), 1e-300)
            push!(rows, (source=s.name, comp=c, Δdeg=cs.Δdeg, specErr=eS, traceErr=eT))
            push!(traces, (source=s.name, comp=c, Δdeg=cs.Δdeg, t=t, cat=uc, ref=ur))
            @printf("  %-24s %s  Δ = %6.2f°   spectrum %.2e   trace (hd %g s) %.2e   (reference disk %.0f s)\n",
                    s.name, c, cs.Δdeg, eS, eT, halfDuration, tref)
        end
    end
    if isdefined(Main, :CairoMakie)
        mkpath(outDir)
        fig = Main.CairoMakie.Figure(size=(1300, 330 * length(sources)))
        for (i, s) in enumerate(sources), (j, c) in enumerate(comps)
            tr = only(filter(x -> x.source == s.name && x.comp == c, traces))
            ax = Main.CairoMakie.Axis(fig[i, j]; xlabel="time (s)", ylabel="displacement (m)",
                title=@sprintf("%s, %s, Δ = %.1f°", s.name, c, tr.Δdeg))
            Main.CairoMakie.lines!(ax, tr.t, tr.ref; color=:black, linewidth=2.0, label="exact source-side")
            Main.CairoMakie.lines!(ax, tr.t, tr.cat; color=:red, linewidth=1.0, label="catalogue")
            (i == 1 && j == 1) && Main.CairoMakie.axislegend(ax; position=:rt)
        end
        f = joinpath(outDir, "check_" * replace(basename(ndFile), ".nd" => "") * ".png")
        Main.CairoMakie.save(f, fig)
        println("figure -> ", f)
    end
    return rows
end

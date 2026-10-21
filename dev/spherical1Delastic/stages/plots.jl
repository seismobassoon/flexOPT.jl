# dev/spherical1Delastic/stages/plots.jl
# Plot/print functions of notebooks/SphericalTI_SH_oct2026.ipynb (cells
# sh-step2-10, sh-step2-ratio, sh-step3-03, sh-step3b-04), unchanged except
# display(fig) -> devSave(fig, ...) (PNG under outputs/<step>/).

function plotSHConvergence(rows; title="")
    fig = Figure(size=(1000, 420))
    ax1 = Axis(fig[1, 1]; xscale=log10, yscale=log10, xlabel="h (km)",
               ylabel="max|A W| / max|A_Ω² W|", title="residual of the exact solution " * title)
    ax2 = Axis(fig[1, 2]; xscale=log10, yscale=log10, xlabel="h (km)",
               ylabel="‖A - sK‖/‖sK‖", title="distance to conventional DSM " * title)
    names = unique(getfield.(rows, :config))
    for nm in names
        rr = filter(r -> r.config == nm, rows)
        h = getfield.(rr, :h_km)
        scatterlines!(ax1, h, getfield.(rr, :resFlexOPT); label=nm)
        scatterlines!(ax2, h, getfield.(rr, :relDiffDSM); label=nm)
    end
    rr = filter(r -> r.config == first(names), rows)
    h = getfield.(rr, :h_km)
    lines!(ax1, h, getfield.(rr, :resDSM); color=:black, linewidth=3, label="DSM conventional")
    h0, r0 = h[end], rr[end].resDSM
    lines!(ax1, h, r0 .* (h ./ h0) .^ 2; color=:gray, linestyle=:dash, label="h²")
    lines!(ax1, h, r0 .* (h ./ h0) .^ 4; color=:gray, linestyle=:dot, label="h⁴")
    axislegend(ax1; position=:rb)
    axislegend(ax2; position=:rb)
    devSave(fig, "convergence " * title)
    return fig
end

# Error constant relative to the conventional FD scheme, at the same h
function shErrorRatios(rows; reference="convFD3")
    for gp in unique(getfield.(rows, :gridParam))
        rr = filter(r -> r.gridParam == gp, rows)
        ref = only(filter(r -> r.config == reference, rr)).resFlexOPT
        @printf("h = %8.3f km   ", rr[1].h_km)
        for r in rr
            r.config == reference && continue
            @printf("%s/%s = %.3e   ", r.config, reference, r.resFlexOPT / ref)
        end
        @printf("DSMconv/%s = %.3e\n", reference, rr[1].resDSM / ref)
    end
end

function plotShellEigen(rows; modes=(1, 5), title="")
    fig = Figure(size=(1100, 450))
    for (ip, im) in enumerate(modes)
        ax = Axis(fig[1, ip]; xscale=log10, yscale=log10, xlabel="h (km)",
                  ylabel="|δΩ/Ω|", title="overtone n=$(im-1) " * title)
        for key in unique([(r.scheme, r.closure) for r in rows])
            rr = filter(r -> (r.scheme, r.closure) == key, rows)
            length(rr[1].relErr) >= im || continue
            h = getfield.(rr, :h_km)
            e = [abs(r.relErr[im]) for r in rr]
            ls = key[1] == "DSMconv" ? :solid : (occursin("traction", key[2]) ? :dash : :solid)
            scatterlines!(ax, h, max.(e, 1e-16); label="$(key[1]) $(key[2])", linestyle=ls)
        end
        ip == length(modes) && Legend(fig[1, ip+1], ax; framevisible=false, labelsize=10)
    end
    devSave(fig, "shell_eigen " * title)
    return fig
end

function plotRecordSection(traces, θdeg; keys=("DSMconv_ref", "DSMconv", "OPT3", "OPT5"), tmax=2000.0, title="")
    fig = Figure(size=(1100, 900))
    ax = Axis(fig[1, 1]; xlabel="time (s)", ylabel="epicentral distance (°)", title="SH (u_φ), f_y = 1 at 590 km, PREM " * title)
    ref = traces[first(keys)]
    sel = ref.t .<= tmax
    scale = 4.0 / maximum(abs.(ref.traces[sel, :]))
    colors = (:black, :orange, :dodgerblue, :crimson)
    for (ik, k) in enumerate(keys)
        tr = traces[k]
        for (ir, θ) in enumerate(θdeg)
            lines!(ax, tr.t[sel], θ .+ scale .* tr.traces[sel, ir]; color=colors[ik], linewidth=ik == 1 ? 2 : 1,
                   label=ir == 1 ? k : nothing)
        end
    end
    axislegend(ax; position=:rb)
    devSave(fig, "record_section " * title); return fig
end

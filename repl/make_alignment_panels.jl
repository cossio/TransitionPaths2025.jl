# Generate the alignment panels of Fig. 2B, Fig. 3B and Fig. 4B.
#
#     julia --project=. repl/make_alignment_panels.jl [outdir]
#
# Writes Fig2B.pdf, Fig3B.pdf, Fig4B.pdf (and .png) into `outdir`
# (default: `figures/`).

import CairoMakie
import TransitionPaths2025 as TP

outdir = isempty(ARGS) ? joinpath(@__DIR__, "..", "figures") : ARGS[1]
mkpath(outdir)

for (name, panel) = [
    ("Fig2B", TP.figure2B_panel),
    ("Fig3B", TP.figure3B_panel),
    ("Fig4B", TP.figure4B_panel),
]
    fig = panel()
    for ext = ("pdf", "png")
        path = joinpath(outdir, "$name.$ext")
        CairoMakie.save(path, fig; px_per_unit = 4)
        println("wrote $path")
    end
end

# Dump the content of the Fig. 2B/3B/4B panels as TSV, to diff them against
# other renderings of the same panels (e.g. the hand-made slides).
#
#     julia --project=. repl/dump_alignment_panels.jl > panels.tsv
#
# One line per row of each panel:
#     panel <TAB> Seq# <TAB> Index <TAB> letter:fill:textcolour × 31
# where `letter` is the residue actually drawn (empty when the panel leaves the
# cell blank), and `fill` / `textcolour` are hex strings.

import TransitionPaths2025 as TP

function dump_panel(name, seq_ids, path)
    seqs = TP.probed_alignment_sequences(seq_ids)
    index = sort(TP.path_index(seqs, path); by = i -> ismissing(i) ? typemax(Int) : i)
    m, n = length(seqs), length(first(seqs))
    chars = [seqs[i][j] for i = 1:m, j = 1:n]
    fills = [TP.alignment_panel_fills(chars[:, j]) for j = 1:n]

    for i = 1:m
        cells = String[]
        for j = 1:n
            drawn = (i == 1 || i == m || chars[i, j] ≠ chars[i - 1, j]) ? string(chars[i, j]) : ""
            fill = fills[j][i]
            text = isempty(drawn) ? "" : uppercase(hex(TP.alignment_panel_text_color(fill)))
            push!(cells, join((drawn, uppercase(hex(fill)), text), ":"))
        end
        println(name, "\t", seq_ids[i], "\t", index[i], "\t", join(cells, "\t"))
    end
end

hex(c) = TP.Makie.Colors.hex(c)

dump_panel("Fig2B", TP.FIG2B_SEQ_IDS, TP.sampled_path_1to1_20240703())
dump_panel("Fig3B", TP.FIG3B_SEQ_IDS, TP.sampled_path_1to4_20240703())
dump_panel("Fig4B", TP.FIG4B_SEQ_IDS, TP.sampled_path_1to2rep1_20240703())

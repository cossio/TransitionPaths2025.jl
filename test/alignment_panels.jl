import TransitionPaths2025
using Test: @test

const TP = TransitionPaths2025

# a column that never changes is fully grey
@test all(==(TP.PANEL_FILL_INVARIANT), TP.alignment_panel_fills(collect("AAAA")))

# endpoints get their own colours, intermediates a distinct one
fills = TP.alignment_panel_fills(collect("ABBC"))
@test fills[1] == TP.PANEL_FILL_FIRST
@test fills[4] == TP.PANEL_FILL_LAST
@test fills[2] == fills[3] == first(TP.PANEL_FILL_INTERMEDIATE)

# a residue going back to the initial one is coloured as the initial one
@test TP.alignment_panel_fills(collect("ABAC"))[3] == TP.PANEL_FILL_FIRST

# two different intermediates get different colours
inter = TP.alignment_panel_fills(collect("ABCD"))
@test inter[2] ≠ inter[3]

# sequences shown in the three panels
for ids = (TP.FIG2B_SEQ_IDS, TP.FIG3B_SEQ_IDS, TP.FIG4B_SEQ_IDS)
    seqs = TP.probed_alignment_sequences(ids)
    @test length(seqs) == length(ids)
    @test all(==(31), length.(seqs))
end

@test only(TP.probed_alignment_sequences(["Seq02"])) == "LPPGWERRADSL-GRTYYVDHNTRTTTWTRP"

# index along the path
@test TP.path_index(["b", "d"], ["a", "b", "c", "d"]) == [2, 4]
@test ismissing(only(TP.path_index(["z"], ["a", "b"])))

# |w|² strips have one value per site
@test length(TP.weight_strip_values(38)) == 31

# the panels build
@test !isnothing(TP.figure2B_panel())
@test !isnothing(TP.figure3B_panel())
@test !isnothing(TP.figure4B_panel())

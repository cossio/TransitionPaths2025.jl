import TransitionPaths2025
using Test: @test

const TP = TransitionPaths2025

# the fills are the ones of the hand-made panels (checked cell by cell against them)
rgb255(c) = round.(Int, 255 .* (Float64(c.r), Float64(c.g), Float64(c.b)))
@test rgb255(TP.PANEL_FILL_INVARIANT) == (217, 217, 217)          # D9D9D9
@test rgb255(TP.PANEL_FILL_FIRST) == (255, 199, 206)              # FFC7CE
@test rgb255(TP.PANEL_FILL_LAST) == (198, 239, 206)               # C6EFCE
@test rgb255(TP.PANEL_FILL_INTERMEDIATE[1]) == (189, 215, 238)    # BDD7EE
@test rgb255(TP.PANEL_FILL_INTERMEDIATE[2]) == (255, 230, 153)    # FFE699

@test rgb255(TP.PANEL_TEXT_FIRST) == (156, 0, 6)                  # 9C0006
@test rgb255(TP.PANEL_TEXT_LAST) == (0, 97, 0)                    # 006100
@test rgb255(TP.PANEL_TEXT_OTHER) == (0, 0, 0)
@test rgb255(TP.PANEL_DOT) == (238, 34, 12)                       # EE220C

# a letter takes the colour of the class of its cell, in every row
@test TP.alignment_panel_text_color(TP.RGBf(TP.PANEL_FILL_FIRST)) == TP.PANEL_TEXT_FIRST
@test TP.alignment_panel_text_color(TP.RGBf(TP.PANEL_FILL_LAST)) == TP.PANEL_TEXT_LAST
@test TP.alignment_panel_text_color(TP.RGBf(TP.PANEL_FILL_INVARIANT)) == TP.PANEL_TEXT_OTHER
@test TP.alignment_panel_text_color(TP.RGBf(TP.PANEL_FILL_INTERMEDIATE[1])) == TP.PANEL_TEXT_OTHER

# a column that never changes is fully grey
@test all(==(TP.PANEL_FILL_INVARIANT), TP.alignment_panel_fills(collect("AAAA")))

# endpoints get their own colours, intermediates a distinct one
fills = TP.alignment_panel_fills(collect("ABBC"))
@test fills[1] == TP.PANEL_FILL_FIRST
@test fills[4] == TP.PANEL_FILL_LAST
@test fills[2] == fills[3] == first(TP.PANEL_FILL_INTERMEDIATE)

# a residue coming back to the initial one has already mutated: it is an intermediate
@test TP.alignment_panel_fills(collect("ABAC"))[3] ∉ (TP.PANEL_FILL_FIRST, TP.PANEL_FILL_LAST)

# two different intermediates get different colours
inter = TP.alignment_panel_fills(collect("ABCD"))
@test inter[2] ≠ inter[3]

# Fig. 4B site 12: the first N is left again for R, so only the second N is final
fills12 = TP.alignment_panel_fills(collect("LNRNN"))
@test fills12[1] == TP.PANEL_FILL_FIRST
@test fills12[2] == first(TP.PANEL_FILL_INTERMEDIATE)
@test fills12[3] == TP.PANEL_FILL_INTERMEDIATE[2]
@test fills12[4] == fills12[5] == TP.PANEL_FILL_LAST

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

# strips are labelled with the paper's numbering of the hidden units
@test TP.paper_hidden_unit_label(38) == "|w₁|²"
@test TP.paper_hidden_unit_label(36) == "|w₂|²"

# the panels build
@test !isnothing(TP.figure2B_panel())
@test !isnothing(TP.figure3B_panel())
@test !isnothing(TP.figure4B_panel())

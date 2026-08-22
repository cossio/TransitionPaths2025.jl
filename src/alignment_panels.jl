###############################################################################
# alignment_panels.jl
#
# Panel-B style alignment tables (Fig. 2B, Fig. 3B, Fig. 4B).
#
# One row per experimentally probed sequence along a path. A residue cell is
# coloured by which endpoint it agrees with (first = red/pink, last = green,
# anything else = a distinct intermediate colour), and letters are only spelled
# out where a mutation happens, except on the two endpoint rows which are shown
# in full. Optionally, strips with the per-site squared norm of the weights of
# the specificity-related hidden units are stacked underneath.
#
# This reproduces the panels that used to be assembled by hand in
# `Seq_notation241211-selected_paths.xlsx` (artifact `2024-12-13_Sequences`).
###############################################################################

# Cell fills. These are the values used in the hand-made spreadsheet version of
# the panels, so that the generated figures match the published ones.
const PANEL_FILL_INVARIANT = colorant"#D9D9D9"   # residue conserved along the path
const PANEL_FILL_FIRST = colorant"#FFC7CE"       # residue as in the first sequence
const PANEL_FILL_LAST = colorant"#C6EFCE"        # residue as in the last sequence
const PANEL_FILL_INTERMEDIATE = (                # neither endpoint, in order of appearance
    colorant"#BDD7ED", colorant"#FFEB84", colorant"#D9C2E9", colorant"#F8CBAD"
)

const PANEL_TEXT_FIRST = colorant"#FF0000"
const PANEL_TEXT_LAST = colorant"#008000"

"""
    alignment_panel_fills(column)

Fill colour of every cell of one alignment column, following the rule used in
Fig. 2B/3B/4B: a column where nothing changes is grey; otherwise a residue
equal to the one of the first (resp. last) sequence takes the first (resp.
last) colour, a residue repeating the one above keeps its colour, and any new
residue takes the next intermediate colour.
"""
function alignment_panel_fills(column::AbstractVector{Char})
    if allequal(column)
        return fill(RGBf(PANEL_FILL_INVARIANT), length(column))
    end

    fills = Vector{RGBf}(undef, length(column))
    intermediates = 0
    for (i, aa) = enumerate(column)
        if aa == first(column)
            fills[i] = PANEL_FILL_FIRST
        elseif aa == last(column)
            fills[i] = PANEL_FILL_LAST
        elseif i > 1 && aa == column[i - 1]
            fills[i] = fills[i - 1]
        else
            intermediates += 1
            fills[i] = PANEL_FILL_INTERMEDIATE[mod1(intermediates, length(PANEL_FILL_INTERMEDIATE))]
        end
    end
    return fills
end

# "38" -> "₃₈", used to label the weight strips
_subscript(n::Integer) = map(c -> Char(0x2080 + (c - '0')), collect(string(n))) |> String

"""
    weight_strip_values(unit; rbm = Eugenio_RBM_20230419(:global))

Per-site squared norm ``|w_i|^2`` of the weights attached to hidden unit `unit`,
i.e. the quantity shown in the strips below the alignment.
"""
function weight_strip_values(unit::Integer; rbm = Eugenio_RBM_20230419(:global))
    return vec(sum(abs2, rbm.w[:, :, unit]; dims=1))
end

"""
    alignment_panel(seqs; labels, index, weight_units, ...)

Build the Fig. 2B / 3B / 4B alignment panel for the aligned sequences `seqs`
(all of the same length, gaps as `-`), ordered along the path.

Keyword arguments:

  * `labels`: row labels, shown in a `Seq#` column (e.g. `["Seq02", ...]`).
  * `index`: position of each sequence along the full sampled path, shown in an
    `Index` column. Use [`path_index`](@ref) to compute it.
  * `weight_units`: hidden units whose `|w|²` profile is drawn underneath, or
    `nothing` for no strips. Defaults to the two specificity-related units.
  * `strip_labels`: labels of those strips; defaults to `|w²₃₈|`, `|w²₃₆|`.
  * `rbm`: RBM providing the weights.
  * `cell`, `fontsize`: size of a residue cell in points, and text size.
"""
function alignment_panel(
    seqs::AbstractVector{<:AbstractString};
    labels::Union{Nothing,AbstractVector} = nothing,
    index::Union{Nothing,AbstractVector} = nothing,
    weight_units = (38, 36),
    strip_labels::Union{Nothing,AbstractVector} = nothing,
    rbm = nothing,
    cell::Real = 20,
    fontsize::Real = 11,
    header_fontsize::Real = 10,
)
    @assert allequal(length, seqs)
    isnothing(labels) || @assert length(labels) == length(seqs)
    isnothing(index) || @assert length(index) == length(seqs)

    m = length(seqs)        # rows (sequences)
    n = length(first(seqs)) # columns (sites)
    chars = [seqs[i][j] for i = 1:m, j = 1:n]

    col_fills = [alignment_panel_fills(chars[:, j]) for j = 1:n]
    fills = [col_fills[j][i] for i = 1:m, j = 1:n]
    variable = [!allequal(chars[:, j]) for j = 1:n]

    units = isnothing(weight_units) ? () : Tuple(weight_units)
    strips = [weight_strip_values(u; rbm = something(rbm, Eugenio_RBM_20230419(:global))) for u = units]
    strip_names = isnothing(strip_labels) ? ["|w²$(_subscript(u))|" for u = units] : strip_labels

    # Layout, in cell units. Row i covers y ∈ [i-1, i]; site j covers x ∈ [j-1, j].
    # The label columns live at negative x, the strips below the last row.
    x_label = isnothing(labels) ? 0.0 : -2.7
    x_index = isnothing(index) ? 0.0 : -0.9
    x_left = min(isnothing(labels) ? 0.0 : x_label - 1.3, isnothing(index) ? 0.0 : x_index - 0.7)
    strip_h = 0.55                    # height of one |w|² strip
    strip_gap = 0.25
    y_strip0 = m + 0.5                # first strip starts here
    y_bottom = isempty(units) ? m : y_strip0 + length(units) * (strip_h + strip_gap)
    x_right = isempty(units) ? n : n + 2.6   # room for the strip labels

    y_top = -0.8                 # room for the site numbers above the first row
    y_end = y_bottom + 0.3

    fig = Makie.Figure()
    ax = Makie.Axis(
        fig[1, 1];
        width = cell * (x_right - x_left),
        height = cell * (y_end - y_top),
        yreversed = true,
    )
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)

    # residue cells
    rects = [Makie.Rect(Makie.Point(j - 1, i - 1), 1, 1) for i = 1:m for j = 1:n]
    Makie.poly!(ax, rects; color = [fills[i, j] for i = 1:m for j = 1:n], strokecolor = :transparent)

    # site numbers on top
    for j = 1:n
        Makie.text!(
            ax, string(j); position = (j - 0.5, -0.35), align = (:center, :center),
            fontsize = header_fontsize
        )
    end

    # residue letters: endpoints in full, intermediates only where they mutate
    for j = 1:n, i = 1:m
        if i == 1
            color = variable[j] ? PANEL_TEXT_FIRST : :black
        elseif i == m
            color = variable[j] ? PANEL_TEXT_LAST : :black
        elseif chars[i, j] == chars[i - 1, j]
            continue
        else
            color = :black
        end
        Makie.text!(
            ax, string(chars[i, j]); position = (j - 0.5, i - 0.5),
            align = (:center, :center), fontsize, color
        )
    end

    # Seq# and Index columns
    if !isnothing(labels)
        Makie.text!(ax, "Seq#"; position = (x_label, -0.35), align = (:center, :center), fontsize = header_fontsize)
        for i = 1:m
            Makie.text!(ax, string(labels[i]); position = (x_label, i - 0.5), align = (:center, :center), fontsize)
        end
    end
    if !isnothing(index)
        Makie.text!(ax, "Index"; position = (x_index, -0.35), align = (:center, :center), fontsize = header_fontsize)
        for i = 1:m
            ismissing(index[i]) && continue   # sequence not on the sampled path
            Makie.text!(ax, string(index[i]); position = (x_index, i - 0.5), align = (:center, :center), fontsize)
        end
    end

    # |w|² strips
    for (s, w2) = enumerate(strips)
        @assert length(w2) == n
        y0 = y_strip0 + (s - 1) * (strip_h + strip_gap)
        lo, hi = extrema(w2)
        colors = [ColorSchemes.get(ColorSchemes.tempo, (w - lo) / (hi - lo)) for w = w2]
        Makie.poly!(
            ax, [Makie.Rect(Makie.Point(j - 1, y0), 1, strip_h) for j = 1:n];
            color = colors, strokecolor = :transparent
        )
        Makie.text!(
            ax, string(strip_names[s]); position = (n + 0.3, y0 + strip_h / 2),
            align = (:left, :center), fontsize
        )
    end

    Makie.xlims!(ax, x_left, x_right)
    Makie.ylims!(ax, y_end, y_top)
    Makie.resize_to_layout!(fig)
    return fig
end

"""
    path_index(seqs, path)

Position of each sequence of `seqs` along `path` (the full sampled path), as
shown in the `Index` column of the panels. `missing` if a sequence is not on
the path.
"""
function path_index(seqs, path)
    idx = indexin(string.(seqs), string.(unique(path)))
    return [isnothing(i) ? missing : i for i = idx]
end

"""
    probed_alignment_sequences(seq_ids)

Aligned sequences of the probed sequences named `seq_ids` (`"Seq02"`, ...), in
the given order. Sequences are looked up in the 2024-12-13 table first, and
then in the 2025-05-31 one (which adds the intermediates of the last batch).
"""
function probed_alignment_sequences(seq_ids::AbstractVector{<:AbstractString})
    df = load_sequences_file_20241213()
    exp = Exp_20240703_sequences()
    # a few sequences are named differently in the two tables
    exp_names = replace(
        exp.names,
        "WW44" => "WW44/147", "WW147" => "WW44/147",
        "WW148" => "WW148/58", "WW149" => "WW9/149"
    )
    aligned = exp.sequences[indexin(df.ww_id, exp_names)]
    by_id = Dict{String,String}(id => string(seq) for (id, seq) = zip(df.seq_id, aligned))

    new = artifact_20250531_new_tested_sequences_load()
    for (id, seq) = zip(new.seq_name, new.sequence)
        get!(by_id, string(id), string(seq))
    end

    unknown = filter(∉(keys(by_id)), seq_ids)
    @assert isempty(unknown) "unknown sequence ids: $unknown"
    return [by_id[id] for id = seq_ids]
end

# Sequences probed along each path, in path order, as shown in the panels.
const FIG2B_SEQ_IDS = ["Seq02", "Seq36", "Seq37", "Seq38", "Seq39", "Seq40", "Seq41", "Seq42"]
const FIG3B_SEQ_IDS = ["Seq02", "Seq49", "Seq50", "Seq51", "Seq47", "Seq52", "Seq53", "Seq54", "Seq35", "Seq48", "Seq27"]
const FIG4B_SEQ_IDS = ["Seq02", "Seq59", "Seq49", "Seq60", "Seq61", "Seq84", "Seq85", "Seq86", "Seq87", "Seq62", "Seq63", "Seq64"]  # Seq84-87 only in the 2025-05-31 batch

"""
    path_alignment_panel(seq_ids, path; sort_by_path = false, kwargs...)

Alignment panel of the probed sequences `seq_ids`, with the `Index` column
giving their position along `path`.

The rows are drawn in the order given by `seq_ids`, which is the order of the
published panels. That is not always the order in which the path visits them:
in Fig. 4B `Seq59` is listed before `Seq49`, while the path goes through
`Seq49` (step 6) before `Seq59` (step 7). Two keywords control what to do
about it:

  * `ascending_index = true` (default) keeps the published row order and lists
    the path positions in increasing order, i.e. `Seq59 → 6`, `Seq49 → 7`, as
    in the panels. All other rows of all three panels are unaffected.
  * `sort_by_path = true` instead reorders the rows by position along the
    path, so that each row differs from the one above by the mutation the path
    actually makes at that step (this puts `Seq49` above `Seq59` in Fig. 4B).
"""
function path_alignment_panel(
    seq_ids, path; sort_by_path::Bool = false, ascending_index::Bool = true, kwargs...
)
    seqs = probed_alignment_sequences(seq_ids)
    index = path_index(seqs, path)
    order = i -> ismissing(i) ? typemax(Int) : i
    if sort_by_path
        perm = sortperm(index; by = order)
        seqs, seq_ids, index = seqs[perm], seq_ids[perm], index[perm]
    elseif ascending_index
        index = sort(index; by = order)
    end
    return alignment_panel(seqs; labels = seq_ids, index, kwargs...)
end

"""
    figure2B_panel(; kwargs...)

Fig. 2B: alignment of the sequences probed along the class I → I path.
"""
figure2B_panel(; kwargs...) =
    path_alignment_panel(FIG2B_SEQ_IDS, sampled_path_1to1_20240703(); kwargs...)

"""
    figure3B_panel(; kwargs...)

Fig. 3B: alignment of the sequences probed along the class I → IV path.
"""
figure3B_panel(; kwargs...) =
    path_alignment_panel(FIG3B_SEQ_IDS, sampled_path_1to4_20240703(); kwargs...)

"""
    figure4B_panel(; kwargs...)

Fig. 4B: alignment of the sequences probed along the class I → II/III path.
This panel carries a single `|w|²` strip, for the hidden unit that separates
class I from class II/III.
"""
figure4B_panel(; weight_units = (38,), kwargs...) =
    path_alignment_panel(FIG4B_SEQ_IDS, sampled_path_1to2rep1_20240703(); weight_units, kwargs...)

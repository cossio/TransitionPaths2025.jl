###############################################################################
# mutational_steps.jl
# Visualise step-wise amino-acid changes in an aligned set of sequences.
###############################################################################

# Build colour lookup per column
function colours_for_column(column::Vector{Char})
	if allequal(column)
		return fill(colorant"#D0D0D0", length(column)) # invariant = grey
	end

	palette = Vector{RGBf}(undef, length(column))
    palette[1] = colorant"#FF8F8F"                     # first, red
    palette[end] = colorant"#A4F2A4"                   # last, green

	tab = ColorSchemes.tab10.colors
	color_index = 0
	for (i, aa) = enumerate(column)
		if aa == first(column)
			palette[i] = colorant"#FF8F8F"  # same as first => red
		elseif aa == last(column)
			palette[i] = colorant"#A4F2A4"   # same as last => green
		elseif aa == column[i-1]
			palette[i] = palette[i-1]
		else
			color_index += 1
			palette[i] = tab[mod1(color_index, length(tab))]
		end
	end

    return palette
end

function aligned_path_plot(seqs::AbstractVector{<:AbstractString})
    @assert allequal(length, seqs)

    m = length(seqs)                       # number of rows
    n = length(seqs[1])                    # number of columns

    chars = [seqs[i][j] for i = 1:m, j = 1:n]   # Char matrix

    col_maps = [colours_for_column(chars[:,j]) for j = 1:n]

    # colour matrix (RGBf)
    cmat = [col_maps[j][i] for i = 1:m, j = 1:n];

    # Per-row mutation counts (relative to previous row)
    diffs = [0; [count(chars[i,k] != chars[i-1,k] for k = 1:n) for i = 2:m]]

    # Plot
    cell = 18                      # pixel size of a square cell
    fig  = Makie.Figure()
    ax   = Makie.Axis(fig[1,1];
        yticks = ((1:m) .- 0.5, string.(1:m)),    # row numbers
        xticks = ((1:n) .- 0.5, string.(1:n)),    # column numbers
        xgridvisible = false, ygridvisible = false,
        yreversed = true,
        xlabel = "", ylabel = "",
		width = cell * n + 60, height= cell * m
    )

    # draw coloured rectangles with residue letters
	xs = repeat(0:n-1, m)
	ys = vcat([fill(i-1, n) for i in 1:m]...)
	rects = Makie.Rect.(Makie.Point.(xs, ys), 1, 1)
	Makie.poly!(ax, rects; color = vec(cmat'), alpha=0.6, strokecolor=:transparent)
	for i = 1:n
		Makie.text!(ax, string(seqs[1][i]); position = (i-0.5, 0.5), align = (:center, :center), fontsize = 12, color=:red, font=:bold)
		Makie.text!(ax, string(seqs[m][i]); position = (i-0.5, m - 0.5), align = (:center, :center), fontsize = 12, color=:green, font=:bold)
		for j = 2:m - 1
			if seqs[j][i] != seqs[j - 1][i]
				Makie.text!(ax, string(seqs[j][i]); position = (i-0.5, j-0.5), align = (:center, :center), fontsize = 12, color=:black)
			end
		end
	end

	Makie.xlims!(ax, 0, n)
	Makie.ylims!(ax, m, 0)

    # right-hand column for mutation counts
    ax2 = Makie.Axis(fig[1,2]; yreversed = true, width=6, height=cell * m)
    for i = 2:m
        Makie.text!(ax2, string(diffs[i]); position = (0.5, i-0.5), align = (:center, :center), fontsize = 13)
    end
	Makie.xlims!(ax2, 0, 1)
	Makie.ylims!(ax2, m, 0)
	Makie.hidedecorations!(ax2)
	Makie.hidespines!(ax2)

	Makie.resize_to_layout!(fig)
    return fig
end







# ---------------------------------------------------------------------------
# 1. Input sequences (all must be the same length, gaps ‘-’ allowed)
seqs = [
    "LPPGWERRADSL-GRTYYVDHNTRTTTWTRP",
    "LPPGWERRVDSL-GRTYYVDHNTRTTTWTRP",
    "LPPGWERRVDPL-GRTYYVDHNTRTTTWTRP",
    "LPPGWERRVDPN-GRTYYVDHNTRTTTWTRP",
    "LPPGWERRVDPN-GRTYYVDHNTRTTTWTRP",
    "LPPGWERRVDPN-GRTYYVDHNTRTTTWTRP",
    "LPPGWERRVDPN-GRVYYVDHNTRTTTWTRP",
    "LPPGWERRVDPN-GRVYYVDHNTRTTTWTDP",
    "LPPGWERRVDPN-GRVYFVDHNTRTTTWTDP",
    "LPPGWERRYDPN-GRVYFVDHNTRTTTWTDP",
    "LPPGWERRYTPN-GRVYFVDHNTRTTTWTDP",
    "LPPGWEIRYTPN-GRVYFVDHNTRTTTWTDP",
    "LPPGWEIRYTPE-GRVYFVDHNTRTTTWTDP",
    "LPPGWEIRYTPE-GRVYFVDHNTRTTTFTDP",
    "LPPGWEIRYTPE-GRRYFVDHNTRTTTFTDP",
    "LPPGWEIRYTPE-GVRYFVDHNTRTTTFTDP",
    "LPPGWEIRYTPE-GVRYFVDHNTRTTTFKDP",
    "LPEGWEIRYTPE-GVRYFVDHNTRTTTFKDP",
    "LPEGWEIRYTRE-GVRYFVDHNTRTTTFKDP"
];

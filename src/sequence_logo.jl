"""
    seqlogo_entropic(p; max_ylim=true)

Plots a conservation sequence logo of protein sequences. Here `p` is a matrix of
probabilities of each letter at each position.
"""
function seqlogo_entropic(p::AbstractMatrix; max_ylim=true)
    @assert size(p, 1) == 21 # nucleotides + gap
    Hmax = log2(size(p, 1))
    w = p ./ sum(p; dims=1)
    H = sum(-xlog2x.(w); dims=1)
    @assert all(0 .≤ H .≤ Hmax)

    #color_scheme = "dmslogo_funcgroup"
    color_scheme = sequence_logo_color_scheme()
    AAs = replace(string(AMINO_ACIDS), '-' => '⊟')

    cons = w .* (Hmax .- H)
    logo = Logomaker.Logo(cons, collect(AAs); color_scheme, flip_below=false)
    max_ylim && logo.ax.set_ylim(0, Hmax)
    logo.ax.set_ylabel("conservation (bits)")
    logo.ax.set_xlabel("site")

    return logo
end

"""
    seqlogo_fields(p; max_ylim=true)

Plots a energetic sequence logo of protein sequences. Here `w` is a matrix of
minus energies of each letter at each position.
"""
function seqlogo_fields(w::AbstractMatrix)
    #color_scheme = "dmslogo_funcgroup"
    color_scheme = sequence_logo_color_scheme()

    AAs = replace(string(AMINO_ACIDS), '-' => '⊟')
    logo = Logomaker.Logo(w, collect(AAs); color_scheme, flip_below=false)
    logo.ax.set_ylabel("fields")
    logo.ax.set_xlabel("site")
    return logo
end

xlog2x(x) = xlogx(x) / log(oftype(x,2))

function sequence_logo_color_scheme()
    Logomaker.color_scheme(
        'C' => "green",
        (a => "orange" for a = ('F', 'W', 'Y'))...,
        (a => "purple" for a = ('Q', 'N', 'S', 'T'))...,
        (a => "black" for a = ('V', 'L', 'I', 'M'))...,
        (a => "blue" for a = ('K', 'R', 'H'))...,
        (a => "red" for a = ('D', 'E'))...,
        (a => "grey" for a = ('A', 'P', 'G'))...,
        '⊟' => "black"
    )
end

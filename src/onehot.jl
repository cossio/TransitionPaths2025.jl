# We encode amino-acid letters in this order (same as Jerome's PGM)
const AMINO_ACIDS = aa"ACDEFGHIKLMNPQRSTVWY-"

function onehot(seq::LongAA)
    seq_ = collect(seq)
    return reshape(seq_, 1, size(seq_)...) .== collect(AMINO_ACIDS)
end

function onehot(seqs::AbstractVector{<:LongAA})
    L = only(unique(length.(seqs)))
    return reshape(mapreduce(onehot, hcat, seqs), 21, L, length(seqs))
    #return reshape(reduce(hcat, onehot.(seqs)), 21, L, length(seqs))
end

onehot(A::AbstractArray{Int8}) = reshape(A, 1, size(A)...) .== 1:21

function aaseq(X::Union{BitMatrix, BitArray{3}})
    @assert size(X, 1) == 21
    return aaseq(potts(X))
end

aaseq(P::AbstractVector{Int8}) = LongAA([AMINO_ACIDS[i] for i in P])
aaseq(P::AbstractMatrix{Int8}) = [aaseq(view(P,:,n)) for n in axes(P, 2)]

potts(X::BitMatrix) = vec(Int8.(first.(Tuple.(argmax(X; dims=1)))))
potts(X::BitArray{3}) = reshape(Int8.(first.(Tuple.(argmax(X; dims=1)))), size(X,2), size(X,3))
potts(s::Union{LongAA, AbstractVector{<:LongAA}}) = potts(onehot(s))

using Test: @test, @testset, @inferred
using TransitionPaths2025: hamming
using TransitionPaths2025: onehot
using BioSequences: @aa_str, LongAA
using Random: randstring

@testset "hamming" begin
    @test hamming(onehot(aa"AGSC"), onehot(aa"AGSC")) == 0
    @test hamming(onehot(aa"AGSC"), onehot(aa"ACSC")) == 1
    @test hamming(onehot(aa"AASC"), onehot(aa"RGSC")) == 2

    seq1 = LongAA(randstring("ACGS", 20))
    seq2 = LongAA(randstring("ACGS", 20))
    @test hamming(seq1, seq2) == hamming(onehot(seq1), onehot(seq2))

    seqs1 = [LongAA(randstring("ACGS", 20)) for _ in 1:100]
    seqs2 = [LongAA(randstring("ACGS", 20)) for _ in 1:100]
    @test hamming(seqs1, seqs2) == hamming(onehot(seqs1), onehot(seqs2))
    @test hamming(seqs1) == hamming(onehot(seqs1))
    @test hamming(seqs1, seqs2[1]) == hamming(onehot(seqs1), onehot(seqs2[1]))
    @test hamming(seqs1[1], seqs2) == hamming(onehot(seqs1[1]), onehot(seqs2))
end

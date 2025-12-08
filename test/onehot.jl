using Test: @test, @testset, @inferred
using BioSequences: @aa_str, LongAA
using LinearAlgebra: I
using TransitionPaths2024: onehot, potts, AMINO_ACIDS

@testset "onehot" begin
    @test onehot(AMINO_ACIDS) == I
    @test potts(AMINO_ACIDS) == 1:21
    @test onehot(Int8.(1:21)) == onehot(AMINO_ACIDS)
end

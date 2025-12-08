using TransitionPaths2025: PF00397_msa_full
using TransitionPaths2025: PF00397_msa_seed
using BioSequences: @aa_str
using BioSequences: LongAA
using Test: @test
using Test: @test_broken
using Test: @testset

@testset "Pfam MSA" begin
    @test all(length.(PF00397_msa_full()) .== 31)
    @test_broken all(length.(PF00397_msa_seed()) .== 31)
end

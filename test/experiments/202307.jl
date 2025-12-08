using Test: @test, @testset
using TransitionPaths2024: Experiments_2307_excel, Exp_202307_results, Probed_Sequences_20231021_path, Exp_202307_sequences

@testset "2307" begin
    @test isfile(Experiments_2307_excel())
    @test size(Exp_202307_results()) == (77, 27)

    @test isfile(Probed_Sequences_20231021_path())
    @test size(Exp_202307_sequences()) == (112, 7)
    @test only(unique(length.(Exp_202307_sequences()[!, "Aligned AA Sequence"]))) == 31
end

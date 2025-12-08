import TransitionPaths2025
using Test: @test

hu_df = TransitionPaths2025.Hu2004_sequences_20230529_aligned()
@test all(==(31), map(length, hu_df.aligned))

df = TransitionPaths2025.Hu2004_sequences_20230529()
aln = TransitionPaths2025.align_to_PF00397(df.WW; remove_inserts=true)
@test aln == hu_df.aligned

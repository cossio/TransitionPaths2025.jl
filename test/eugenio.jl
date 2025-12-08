import FASTX
import HDF5
using RestrictedBoltzmannMachines: energy
using RestrictedBoltzmannMachines: free_energy
using Test: @test, @testset
using BioSequences: LongAA
using TransitionPaths2024: Eugenio_fasta_20230419_path
using TransitionPaths2024: Eugenio_MSA_20230419
using TransitionPaths2024: Eugenio_Probed_Sequences_202306
using TransitionPaths2024: Eugenio_Probed_Sequences_202306_path
using TransitionPaths2024: Eugenio_RBM_20230419
using TransitionPaths2024: Eugenio_RBM_20230419
using TransitionPaths2024: Eugenio_RBM_20230419_loglikelihood
using TransitionPaths2024: Eugenio_RBM_20230419_path
using TransitionPaths2024: Eugenio_RBM_20230419_probed_sequences_free_energies_eval_from_python
using TransitionPaths2024: Eugenio_RBM_logZ_20230419
using TransitionPaths2024: Eugenio_RBM_specific_20230424
using TransitionPaths2024: Eugenio_RBM_typeI_20230424_path
using TransitionPaths2024: Eugenio_RBM_typeII_20230424_path
using TransitionPaths2024: Eugenio_RBM_typeIV_20230424_path
using TransitionPaths2024: Exp_20240703_sequences
using TransitionPaths2024: Literature_sequences_20230424
using TransitionPaths2024: onehot

@test isfile(Eugenio_RBM_20230419_path())
@test isfile(Eugenio_fasta_20230419_path())

@testset "Eugenio's RBM (20230419)" begin
    rbm = Eugenio_RBM_20230419()
    v = HDF5.h5open(Eugenio_RBM_20230419_path(), "r") do hdf
        test_data_v = HDF5.read(hdf, "test_data_v")
        test_data_h = HDF5.read(hdf, "test_data_h")
        test_data_f = HDF5.read(hdf, "test_data_f")
        test_data_E = HDF5.read(hdf, "test_data_E")

        # onehot encode
        v = reshape(test_data_v, 1, size(test_data_v)...) .== 0:20

        @test energy(rbm, v, test_data_h) ≈ test_data_E rtol=1e-5
        @test free_energy(rbm, v) ≈ test_data_f rtol=1e-3
    end
end

@testset "Eugenio's MSA (20230419)" begin
    msa = Eugenio_MSA_20230419()
    @test all(length.(msa) .== 31)
end

@testset "Eugenio's literature sequences (20230424)" begin
    data = Literature_sequences_20230424()
    @test all(length.(data.sequences_I) .== 31)
    @test all(length.(data.sequences_II) .== 31)
    @test all(length.(data.sequences_III) .== 31)
    @test all(length.(data.sequences_IV) .== 31)
end

@testset "Eugenio probed sequences (202306)" begin
    @test isfile(Eugenio_Probed_Sequences_202306_path())
    seqs = Eugenio_Probed_Sequences_202306()
    @test length(seqs) == 32
    @test only(unique(length.(seqs))) == 31
end

@testset "Eugenio RBMs scores" begin
    sequences = onehot(LongAA.([
        "LPPGWEKRMSRS-GRVYYVNHITRASQWERP",
        "LPPGWEKRMSRS-GRVYYVNHITNASQWERP",
        "LPPGWEKRMSRSSGRVYYVNHITRASQWERP",
        "LPPGWEKRMSRSSGRVYYVNHITNASQWERP",
    ]))

    @test Eugenio_RBM_20230419_loglikelihood(:I, sequences) ≈ [-34.07231, -39.860474, -37.56433, -43.245316] rtol=0.005
    @test Eugenio_RBM_20230419_loglikelihood(:II, sequences) ≈ [-42.608047, -42.158203, -45.2482, -44.420837] rtol=0.005
    @test Eugenio_RBM_20230419_loglikelihood(:IV, sequences) ≈ [-23.032043, -18.61139, -16.968628, -11.998077] rtol=0.005
    @test Eugenio_RBM_20230419_loglikelihood(:global, sequences) ≈ [-27.112946,-29.897995, -28.901306, -30.94751] rtol=0.005
end

@testset "Eugenio RBMs Eeff" begin
    @test free_energy(Eugenio_RBM_20230419(:I), onehot(Exp_20240703_sequences().sequences)) ≈ Eugenio_RBM_20230419_probed_sequences_free_energies_eval_from_python(:I)
    @test free_energy(Eugenio_RBM_20230419(:II), onehot(Exp_20240703_sequences().sequences)) ≈ Eugenio_RBM_20230419_probed_sequences_free_energies_eval_from_python(:II)
    @test free_energy(Eugenio_RBM_20230419(:IV), onehot(Exp_20240703_sequences().sequences)) ≈ Eugenio_RBM_20230419_probed_sequences_free_energies_eval_from_python(:IV)
    @test free_energy(Eugenio_RBM_20230419(:global), onehot(Exp_20240703_sequences().sequences)) ≈ Eugenio_RBM_20230419_probed_sequences_free_energies_eval_from_python(:global)
end

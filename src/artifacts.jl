# PF00397 alignment, downloaded from Pfam on 20230419
PF00397_seed_20230419_path() = joinpath(artifact"20230419_PF00397", "PF00397.alignment.seed")
PF00397_full_20230419_path() = joinpath(artifact"20230419_PF00397", "PF00397.alignment.full")

# Eugenio's data, 20230419
Eugenio_RBM_20230419_path() = joinpath(artifact"20230419_Eugenio", "WW_RBM.hdf5")
Eugenio_fasta_20230419_path() = joinpath(artifact"20230419_Eugenio", "WW_domain_MSA.fasta")

# Eugenio's specific RBMs, 20230424
Eugenio_RBM_typeI_20230424_path() = joinpath(artifact"20230424_Eugenio_Specific_RBMs", "WW_RBM_type1.hdf5")
Eugenio_RBM_typeII_20230424_path() = joinpath(artifact"20230424_Eugenio_Specific_RBMs", "WW_RBM_type2.hdf5")
Eugenio_RBM_typeIV_20230424_path() = joinpath(artifact"20230424_Eugenio_Specific_RBMs", "WW_RBM_type4.hdf5")

# Sequnces validated in the literature, originally collected by Jerome, 20230424
Literature_sequences_20230424_path() = joinpath(artifact"Eugenio_20230424_Tested_sequences_from_literature", "WW_test_sequences.hdf5")

# Labeled sequences from Hu 2004 paper, "A map of WW domain family interactions"
Hu2004_sequences_20230529_path() = joinpath(artifact"20230425_WW_Hu_2004", "WW_2004.tsv")

# Labeled sequences from Hu 2004 paper, "A map of WW domain family interactions"
Eugenio_Probed_Sequences_202306_path() = joinpath(artifact"Eugenio_Probed_Sequences_202306", "32_for_marco.faa")

# Experiments including AdvRBM sequences, from 2023-09.
Experiments_2307_excel() = joinpath(artifact"Experiments_2307", "analysis 2307.xlsx")

# Probed sequences in experimenst, up to 2023-10-21
Probed_Sequences_20231021_path() = joinpath(artifact"Probed_Sequences_20231021", "WW-DATA - Fixed.tsv")

# Probed sequences, with alignments, up to 2024-07-03
Aligned_Probed_Sequences_20240703_path() = joinpath(artifact"Aligned_Probed_Sequences_20240703", "WW-DATA_aln_with_rev.xlsx")

Probed_Sequences_Eeff_from_Python_20240703_path() = joinpath(artifact"Probed_Sequences_Eeff_from_Python_20240703")

# Repository with Eugenio's work (paths, RBMs, etc.). Original source: https://github.com/eugeniomauri1/trans_path
Transition_Paths_Eugenio_Repo_20240704_path() = artifact"Transition_Paths_Eugenio_Repo_20240704"

function artifact_20241213_sequences_file()
    return joinpath(artifact"2024-12-13_Sequences", "Seq_notation241211-selected_paths.xlsx")
end

function artifact_20250531_new_tested_sequences_file()
    return joinpath(artifact"2025-05-31_New_tested_sequences", "AR250527.xlsx")
end

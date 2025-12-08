function sampled_path_1to1_20240703()
    fasta_file = joinpath(Transition_Paths_Eugenio_Repo_20240704_path(), "full_paths", "1to1_batch2_full.faa")
    return LongAA.(FASTX.sequence.(FASTX.FASTA.Reader(open(fasta_file))))
end

function sampled_path_1to2rep2_20240703()
    fasta_file = joinpath(Transition_Paths_Eugenio_Repo_20240704_path(), "full_paths", "1to2(rep2)_batch4_full.faa")
    return LongAA.(FASTX.sequence.(FASTX.FASTA.Reader(open(fasta_file))))
end

function sampled_path_1to2rep1_20240703()
    fasta_file = joinpath(Transition_Paths_Eugenio_Repo_20240704_path(), "full_paths", "1to2(rep1)_batch4_full.faa")
    return LongAA.(FASTX.sequence.(FASTX.FASTA.Reader(open(fasta_file))))
end

function sampled_path_1to4_20240703()
    fasta_file = joinpath(Transition_Paths_Eugenio_Repo_20240704_path(), "full_paths", "1to4_batch3_full.faa")
    return LongAA.(FASTX.sequence.(FASTX.FASTA.Reader(open(fasta_file))))
end

function sampled_path_1to4direct_20240703()
    fasta_file = joinpath(Transition_Paths_Eugenio_Repo_20240704_path(), "full_paths", "1to4dir_batch3_full.faa")
    return LongAA.(FASTX.sequence.(FASTX.FASTA.Reader(open(fasta_file))))
end

function sampled_path_1to4batch1_20240703()
    fasta_file = joinpath(Transition_Paths_Eugenio_Repo_20240704_path(), "full_paths", "1to4_batch1_full.faa")
    return LongAA.(FASTX.sequence.(FASTX.FASTA.Reader(open(fasta_file))))
end

function sampled_path_1to4batch2_20240703()
    fasta_file = joinpath(Transition_Paths_Eugenio_Repo_20240704_path(), "full_paths", "1to4_batch2_full.faa")
    return LongAA.(FASTX.sequence.(FASTX.FASTA.Reader(open(fasta_file))))
end

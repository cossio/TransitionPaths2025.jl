function PF00397_msa_full()
    path = PF00397_full_20230419_path()
    fasta = HMMER.esl_reformat("afa", path; informat="stockholm").o
    seqs = FASTX.sequence.(FASTX.FASTA.Reader(open(fasta)))
    seqs = filter.(c -> !(c == '.' || islowercase(c)), seqs) # remove inserts
    return LongAA.(seqs)
end

function PF00397_msa_seed()
    path = PF00397_seed_20230419_path()
    fasta = HMMER.esl_reformat("afa", path; informat="stockholm").o
    seqs = FASTX.sequence.(FASTX.FASTA.Reader(open(fasta)))
    seqs = filter.(c -> !(c == '.' || islowercase(c)), seqs) # remove inserts
    @warn "inserts are annotated differently in the seed, so this alignment is broken"
    return LongAA.(seqs)
end

function PF00397_path(which::Symbol = :full)
    if which === :full
        return PF00397_full_20230419_path()
    elseif which === :seed
        return PF00397_seed_20230419_path()
    else
        throw(ArgumentError("Invalid argument: $which"))
    end
end

# Labeled sequences from Hu 2004 paper, "A map of WW domain family interactions" (not aligned)
Hu2004_sequences_20230529() = DataFrame(CSV.File(Hu2004_sequences_20230529_path(), delim='\t'))

function Hu2004_sequences_20230529_aligned()
    hu2004_df = Hu2004_sequences_20230529()

    # prepare FASTA
    hu2004_fasta = tempname()
    FASTX.FASTA.Writer(open(hu2004_fasta, "w")) do writer
        for (i, seq) in enumerate(hu2004_df.WW)
            write(writer, FASTX.FASTA.Record(string(i), seq))
        end
    end

    # align
    pfam_hmm = Pfam.Pfam_A_hmm()
    hmm = HMMER.hmmfetch(pfam_hmm, "PF00397.29").o
    aln = HMMER.hmmalign(hmm, hu2004_fasta; outformat="afa")

    # load aligned sequences from aln.o output file
    aligned_sequences = String[]
    FASTX.FASTA.Reader(open(aln.o)) do reader
        for record in reader
            push!(aligned_sequences, filter(c -> !islowercase(c) && c != '.', FASTX.sequence(record)))
        end
    end

    # add aligned sequences column to DataFrame
    hu2004_df.aligned = aligned_sequences
    return hu2004_df
end

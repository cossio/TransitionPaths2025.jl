function align_to_PF00397(sequences::AbstractVector; remove_gaps::Bool = false, remove_inserts::Bool)
    # prepare FASTA
    fasta = tempname()
    FASTX.FASTA.Writer(open(fasta, "w")) do writer
        for (i, seq) in enumerate(sequences)
            if remove_gaps
                seq_to_write = replace(seq, '-' => "")
            else
                seq_to_write = seq
            end
            write(writer, FASTX.FASTA.Record(string(i), seq_to_write))
        end
    end

    # align
    pfam_hmm = Pfam.Pfam_A_hmm()
    hmm = HMMER.hmmfetch(pfam_hmm, "PF00397.29").o
    aln = HMMER.hmmalign(hmm, fasta; outformat="afa")

    # load aligned sequences from aln.o output file
    aligned_sequences = String[]
    FASTX.FASTA.Reader(open(aln.o)) do reader
        for record = reader
            if remove_inserts
                push!(aligned_sequences, filter(c -> !islowercase(c) && c != '.', FASTX.sequence(record)))
            else
                push!(aligned_sequences, FASTX.sequence(record))
            end
        end
    end

    return aligned_sequences
end

function Exp_202307_results()
    xls = readxlsx(Experiments_2307_excel())["test 2307"]

    column_labels = vec(xls[1,begin:end - 12])
    column_labels[6:9] .*= "_1"
    column_labels[10:13] .*= "_2"
    column_labels[14:16] .*= "_diff_1"
    column_labels[17:19] .*= "_diff_2"
    column_labels[20:22] .*= "_rel_1"
    column_labels[23:25] .*= "_rel_2"

    df = DataFrame(xls[2:end, begin:end - 12], string.(column_labels))

    # Add column with sequences
    seqs_df = Exp_202307_sequences()
    df.AASequence = map(Returns(""), df[!,"WW#"])
    df.AlignedAASequence = map(Returns(""), df[!,"WW#"])
    for (k, i) = enumerate(indexin(df[!, "WW#"], seqs_df.NAME))
        if !isnothing(i)
            df.AASequence[k] = seqs_df[i, "AA Sequence"]
            df.AlignedAASequence[k] = seqs_df[i, "Aligned AA Sequence"]
        end
    end

    return df
end

function Exp_202307_sequences()
    return CSV.read(Probed_Sequences_20231021_path(), DataFrame)
end

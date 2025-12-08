function load_sequences_file_20241213()
    xls = readxlsx(artifact_20241213_sequences_file())["Selected_paths"]

    rows = [3:29; 32:41; 44:61; 66:85; 90:105; 110:125]

    df = DataFrame(
        ww_id = string.(xls[i, 1] for i = rows),
        seq_id = string.(xls[i, 2] for i = rows),
        notation = string.(xls[i, 3] for i = rows),
        sequence = string.(xls[i, 4] for i = rows),
        C1 = float.(xls[i, 5] for i = rows),
        C2 = float.(xls[i, 6] for i = rows),
        C4 = float.(replace([xls[i, 7] for i = rows], "-" => missing)),
        eC1 = float.(replace([xls[i, 8] for i = rows], "-" => missing)),
        eC2 = float.(replace([xls[i, 9] for i = rows], "-" => missing)),
        eC4 = float.(replace([xls[i, 10] for i = rows], "-" => missing)),
        group = [
            fill("Wild-types", length(3:29));
            fill("1st_batch", length(32:41));
            fill("2nd_batch", length(44:61));
            fill("3rd_batch", length(66:85));
            fill("4th_batch", length(90:105));
            fill("Scrambled_paths", length(110:125))
        ]
    )

    return df
end

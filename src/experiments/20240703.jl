function Exp_20240703_sequences()
    xls = readxlsx(Aligned_Probed_Sequences_20240703_path())["Foglio1"][:]
    names = map(identity, xls[2:151,1]) # names, e.g. WW1, WW2, ...
    sequences = map(LongAA, xls[2:151,5]) # aligned sequences
    return (; names, sequences)
end

function Exp_20240703_sequences_2()
    #xls = readxlsx(Aligned_Probed_Sequences_20240703_path())["Foglio1"][:]
    xls = readxlsx(pkgdir(@__MODULE__, "data", "WW-DATA_aln_with_rev.xlsx"))["Foglio1"][:]
    names = map(identity, xls[2:151,1]) # names, e.g. WW1, WW2, ...
    sequences = map(LongAA, xls[2:151,5]) # aligned sequences
    return (; names, sequences)
end

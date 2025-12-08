function artifact_20250531_new_tested_sequences_load()
    row_range = 2:13
    xls = readxlsx(artifact_20250531_new_tested_sequences_file())
    seq_name = map(identity, xls["Sheet1"][:][row_range, 1])
    notation = map(identity, xls["Sheet1"][:][row_range, 2])
    sequence = map(identity, xls["Sheet1"][:][row_range, 3])
    C1_response = map(identity, xls["Sheet1"][:][row_range, 4])
    C2_response = map(identity, xls["Sheet1"][:][row_range, 5])
    C4_response = map(identity, xls["Sheet1"][:][row_range, 6])
    C1_response_err = map(identity, xls["Sheet1"][:][row_range, 7])
    C2_response_err = map(identity, xls["Sheet1"][:][row_range, 8])
    C4_response_err = map(identity, xls["Sheet1"][:][row_range, 9])
    return (; seq_name, notation, sequence, C1_response, C2_response, C4_response, C1_response_err, C2_response_err, C4_response_err)
end

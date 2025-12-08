function response_data_xls_20240703()
    return pkgdir(@__MODULE__, "data", "2024-07-03-analysis.xlsx")
end

function response_data_xls_20240716()
    return pkgdir(@__MODULE__, "data", "2024-07-16-analysis.xlsx")
end

function response_data_path_1to1_20240703()
    xls = readxlsx(response_data_xls_20240716())
    ww_names = map(identity, xls["Paths"][:][2:9, 1])
    C1_response = map(identity, xls["Paths"][:][2:9, 2])
    C1_response_err = map(identity, xls["Paths"][:][2:9, 3])
    return (; ww_names, C1_response, C1_response_err)
end

function response_data_path_1to2rep1_20240703()
    xls = readxlsx(response_data_xls_20240716())
    ww_names = map(identity, xls["Paths"][:][14:21, 1])
    C1_response = map(identity, xls["Paths"][:][14:21, 2])
    C2_response = map(identity, xls["Paths"][:][14:21, 3])
    C1_response_err = map(identity, xls["Paths"][:][14:21, 4])
    C2_response_err = map(identity, xls["Paths"][:][14:21, 5])
    return (; ww_names, C1_response, C2_response, C1_response_err, C2_response_err)
end

function response_data_path_1to2rep2_20240703()
    xls = readxlsx(response_data_xls_20240716())
    ww_names = map(identity, xls["Paths"][:][14:21, 7])
    C1_response = map(identity, xls["Paths"][:][14:21, 8])
    C2_response = map(identity, xls["Paths"][:][14:21, 9])
    C1_response_err = map(identity, xls["Paths"][:][14:21, 10])
    C2_response_err = map(identity, xls["Paths"][:][14:21, 11])
    return (; ww_names, C1_response, C2_response, C1_response_err, C2_response_err)
end

function response_data_path_1to4_20240703()
    # xls = readxlsx(response_data_xls_20240703())
    # ww_names = map(identity, xls["Paths"][:][26:33, 1])
    # C1_response = map(identity, xls["Paths"][:][26:33, 2])
    # C4_response = map(identity, xls["Paths"][:][26:33, 3])
    # C1_response_err = map(identity, xls["Paths"][:][26:33, 4])
    # C4_response_err = map(identity, xls["Paths"][:][26:33, 5])

    csv = CSV.read(joinpath(artifact"Response_Data_1to4_20240708", "CSV", "1to4-1to4.csv"), DataFrame)

    # WW44 was a duplicate, so we remove it
    ww_names = csv.Name[csv.Name .!= "WW44"]
    C1_response = csv.C1[csv.Name .!= "WW44"]
    C4_response = csv.C4[csv.Name .!= "WW44"]
    C1_response_err = csv.eC1[csv.Name .!= "WW44"]
    C4_response_err = csv.eC4[csv.Name .!= "WW44"]

    return (; ww_names, C1_response, C4_response, C1_response_err, C4_response_err)
end

function response_data_path_1to4direct_20240703()
    # xls = readxlsx(response_data_xls_20240703())
    # ww_names = map(identity, xls["Paths"][:][26:33, 7])
    # C1_response = map(identity, xls["Paths"][:][26:33, 8])
    # C4_response = map(identity, xls["Paths"][:][26:33, 9])
    # C1_response_err = map(identity, xls["Paths"][:][26:33, 10])
    # C4_response_err = map(identity, xls["Paths"][:][26:33, 11])
    # return (; ww_names, C1_response, C4_response, C1_response_err, C4_response_err)

    csv = CSV.read(joinpath(artifact"Response_Data_1to4_20240708", "CSV", "1to4direct-1to4direct.csv"), DataFrame)

    # WW58 and WW09 are duplicates, so we remove them
    ww_names = csv.Name[(csv.Name .!= "WW58") .& (csv.Name .!= "WW09")]
    C1_response = csv.C1[(csv.Name .!= "WW58") .& (csv.Name .!= "WW09")]
    C4_response = csv.C4[(csv.Name .!= "WW58") .& (csv.Name .!= "WW09")]
    C1_response_err = csv.eC1[(csv.Name .!= "WW58") .& (csv.Name .!= "WW09")]
    C4_response_err = csv.eC4[(csv.Name .!= "WW58") .& (csv.Name .!= "WW09")]

    return (; ww_names, C1_response, C4_response, C1_response_err, C4_response_err)
end

function response_data_all_20240703()
    xls = readxlsx(response_data_xls_20240703())

    response_data_C1 = map(identity, xls["compress"][:][4:67, 2])
    response_data_C2 = map(identity, xls["compress"][:][4:67, 3])
    response_data_C4 = map(identity, xls["compress"][:][4:67, 5])
    response_data_C1_err = map(identity, xls["compress"][:][4:67, 8])
    response_data_C2_err = map(identity, xls["compress"][:][4:67, 9])
    response_data_C4_err = map(identity, xls["compress"][:][4:67, 11])

    response_data_WW_names = map(xls["compress"][:][4:67, 1]) do str
        # We correct names like WW07 to WW7
        if str ∈ ["WW0$i" for i = 1:9]
            return replace(str, '0' => "")
        else
            return str
        end
    end

    return (;
        response_data_C1, response_data_C2, response_data_C4,
        response_data_C1_err, response_data_C2_err, response_data_C4_err,
        response_data_WW_names
    )
end

function response_data_path_1to2rev_20240703()
    xls = readxlsx(response_data_xls_20240703())
    ww_names = map(identity, xls["Paths"][:][51:58, 1])
    C1_response = map(identity, xls["Paths"][:][51:58, 2])
    C2_response = map(identity, xls["Paths"][:][51:58, 3])
    C1_response_err = map(identity, xls["Paths"][:][51:58, 4])
    C2_response_err = map(identity, xls["Paths"][:][51:58, 5])
    return (; ww_names, C1_response, C2_response, C1_response_err, C2_response_err)
end

function response_data_path_1to4rev_20240703()
    xls = readxlsx(response_data_xls_20240703())
    ww_names = map(identity, xls["Paths"][:][51:58, 7])
    C1_response = map(identity, xls["Paths"][:][51:58, 8])
    C4_response = map(identity, xls["Paths"][:][51:58, 9])
    C1_response_err = map(identity, xls["Paths"][:][51:58, 10])
    C4_response_err = map(identity, xls["Paths"][:][51:58, 11])
    return (; ww_names, C1_response, C4_response, C1_response_err, C4_response_err)
end

function response_data_all_v2_20240708()
    df = CSV.read(joinpath(artifact"Response_Data_All_20240708", "All-All.csv"), DataFrame)
    name = [replace(p, "WW0" => "WW") for p = df.Peptide]
    return (; name, df.C1, df.C2, C4=df.var"C4.2", df.C1err, df.C2err, C4err=df.var"C4.2err", type=df.Type)
end

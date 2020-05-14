import CSV

function get_detailled_crs_data(path::String, filebase::String, selected_id::Integer)

    # select all relevant files
    files = path .* filter!(s->occursin(r"^" * "$filebase", s), readdir(path))

    # read files
    df = mapreduce(CSV.read, vcat, files)

    # delete all empty rows
    CSV.dropmissing!(df)

    # filter dataframe by relevant id
    df = df[df.ID .== selected_id,:]

    return sort!(df[df[:, :ID] .== selected_id,:])
end

function get_detailled_shock_data(path::String, selected_id::Integer)
    return get_detailled_crs_data(path, "CRs_details_shock_injection", mpi_ranks, selected_id)
end

function get_detailled_Dpp_data(path::String, selected_id::Integer)
    return get_detailled_crs_data(path, "CRs_details_Dpp", mpi_ranks, selected_id)
end

function get_detailled_radiative_data(path::String, selected_id::Integer)
    return get_detailled_crs_data(path, "CRs_details_radiative", mpi_ranks, selected_id)
end

function get_detailled_adiabatic_data(path::String, selected_id::Integer)
    return get_detailled_crs_data(path, "CRs_details_adiabatic", mpi_ranks, selected_id)
end

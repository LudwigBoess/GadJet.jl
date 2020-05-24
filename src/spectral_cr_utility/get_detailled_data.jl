import CSV

function get_detailled_crs_data(path::String, filebase::String)

    # select all relevant files
    files = path .* filter!(s->occursin(r"^" * "$filebase", s), readdir(path))

    # read files
    df = mapreduce(CSV.read, vcat, files)

    # delete all empty rows
    CSV.dropmissing!(df)

    return sort!(df)
end

function get_detailled_shock_data(path::String, selected_id::Integer=-1)
    if (selected_id == -1)
        return get_detailled_crs_data(path, "CRs_details_shock_injection")
    else
        df = get_detailled_crs_data(path, "CRs_details_shock_injection")
        return sort!(df[df[:, :ID] .== selected_id,:])
    end
end

function get_detailled_Dpp_data(path::String, selected_id::Integer=-1)
    if (selected_id == -1)
        return get_detailled_crs_data(path, "CRs_details_Dpp")
    else
        df = get_detailled_crs_data(path, "CRs_details_Dpp")
        return sort!(df[df[:, :ID] .== selected_id,:])
    end
end

function get_detailled_radiative_data(path::String, selected_id::Integer=-1)
    if (selected_id == -1)
        return get_detailled_crs_data(path, "CRs_details_radiative")
    else
        df = get_detailled_crs_data(path, "CRs_details_radiative")
        return sort!(df[df[:, :ID] .== selected_id,:])
    end
end

function get_detailled_adiabatic_data(path::String, selected_id::Integer=-1)
    if (selected_id == -1)
        return get_detailled_crs_data(path, "CRs_details_adiabatic")
    else
        df = get_detailled_crs_data(path, "CRs_details_adiabatic")
        return sort!(df[df[:, :ID] .== selected_id,:])
    end
end

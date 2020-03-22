import CSV.read

function get_detailled_crs_data(path::String, filebase::String, mpi_ranks::Int, selected_id::Int)

    read_first = false

    @inbounds for mpi_rank = 0:(mpi_ranks-1)

        fi = path * filebase * "_$mpi_rank.txt"

        d = CSV.read(fi)

        df2 = d[d.ID .== selected_id,:]

        if length(df2[!, 1]) > 0

            if (read_first == false)
                df = df2
                read_first = true
            else
                append!(df, df2)
            end
        end

    end

    return sort!(df)
end

function get_detailled_shock_data(path::String, mpi_ranks::Int, selected_id::Int)
    return get_detailled_crs_data(path, "CRs_details_shock_injection", mpi_ranks, selected_id)
end

function get_detailled_Dpp_data(path::String, mpi_ranks::Int, selected_id::Int)
    return get_detailled_crs_data(path, "CRs_details_Dpp", mpi_ranks, selected_id)
end

function get_detailled_radiative_data(path::String, mpi_ranks::Int, selected_id::Int)
    return get_detailled_crs_data(path, "CRs_details_radiative", mpi_ranks, selected_id)
end

function get_detailled_adiabatic_data(path::String, mpi_ranks::Int, selected_id::Int)
    return get_detailled_crs_data(path, "CRs_details_adiabatic", mpi_ranks, selected_id)
end

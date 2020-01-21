using CSV

function get_detailled_shock_data(path::String, mpi_ranks::Int, selected_id::Int)

    fi = path * "CRs_details_shock_injection_0.txt"

    d = CSV.read(fi)

    df = d[d.ID .== selected_id,:]

    if mpi_ranks > 1

        for mpi_rank = 1:(mpi_ranks-1)

            fi = path * "CRs_details_shock_injection_$mpi_rank.txt"

            d = CSV.read(fi)

            df2 = d[d.ID .== selected_id,:]

            append!(df, df2)

        end

    end

    return sort!(df)

end

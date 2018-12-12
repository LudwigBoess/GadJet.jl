# ------------------------------------------------------------------------------
# ---------------------- Dictionary functions ----------------------------------
# ------------------------------------------------------------------------------

function head_to_dict(filename::String, cosmic_rays::Bool=false,
                      n_spec_indices::Int64=0)

        header = Dict()

        f = open(filename)
        blocksize = read(f, Int32)

        if blocksize[1] == 8
            swap = 0
            snap_format = 2
        elseif blocksize[1] == 256
            swap = 0
            snap_format = 1
        else
            blocksize[1] = bswap(blocksize[1])
            if blocksize[1] == 8
                swap = 1
                snap_format = 2
            elseif blocksize[1] == 256
                swap = 1
                snap_format = 1
            else
                println("incorrect file format encountered when reading header of", filename)
            end
        end

        if snap_format == 2
            seek(f, 16)
            skip_line = read(f, Int32)
        end

        header["snap_format"] = snap_format
        header["PartTypes"] = ["PartType0", "PartType1", "PartType2",
                               "PartType3", "PartType4", "PartType5"]
        header["npart"] = read!(f, Array{Int32,1}(undef,6))
        header["massarr"] = read!(f, Array{Float64,1}(undef,6))
        header["time"] = read(f, Float64)
        header["redshift"] = read(f, Float64)
        header["flag_sfr"] = read(f, Int32)
        header["flag_feedback"] = read(f, Int32)
        header["nall"] = read!(f, Array{UInt32,1}(undef,6))
        header["flag_cooling"] = read(f, Int32)
        header["num_files"] = read(f, Int32)
        header["boxsize"] = read(f, Float64)
        header["omega_m"] = read(f, Float64)
        header["omega_l"] = read(f, Float64)
        header["hubble"] = read(f, Float64)
        header["flag_stellarage"] = read(f, Int32)
        header["flag_metals"] = read(f, Int32)
        header["npartTotalHighWord"] = read!(f, Array{UInt32,1}(undef,6))
        header["flag_entropy_instead_u"] = read(f, Int32)
        header["flag_doubleprecision"] = read(f, Int32)
        header["flag_ic_info"] = read(f, Int32)
        header["lpt_scalingfactor"] = read(f, Float32)

        # skipline = 0
        # iterations = 1
        # while skipline != Int32(256)
        #    skipline = read(f, Int32)
        #    iterations += 1
        # end
        # println(iterations)
        close(f)

        return header

end


function snap_to_dict(filename::String, try_info::Bool=true)

        data = Dict()

        data["Header"] = head_to_dict(filename)

        if data["Header"]["snap_format"] == 2

            if try_info

                    info = get_info(filename)

                    if typeof(info) == Array{Info_Line,1}

                        data = snap_2_d_info(filename, data, info)

                    else

                        data = snap_2_d(filename, data)

                    end

            else

                data = snap_2_d(filename, data)

            end
        else

            data = snap_1_d_old(filename, data)

        end

        return data

end


function snap_2_d_info(filename::String, d::Dict{Any,Any}, info::Array{Info_Line,1})

    f = open(filename)

    seek(f, 300)

    for i ∈ 1:length(d["Header"]["npart"])
        if d["Header"]["npart"][i] != Int32(0)
            d[d["Header"]["PartTypes"][i]] = Dict()
            d[d["Header"]["PartTypes"][i]]["MASS"] = Float32.(d["Header"]["massarr"][i] .* ones(d["Header"]["npart"][i]))
        end
    end

    for i ∈ 1:length(info)
        #println("Reading block: ", info[i].block_name)
        d = read_block_with_info(f, d, info[i])
        p = position(f)
        seek(f, p+24)
    end

    close(f)

    return d

end

function snap_2_d(filename::String, data::Dict{Any,Any})

    f = open(filename)

    seek(f, 296)

    N = sum(data["Header"]["npart"])
    skipsize = read(f, Int32)
    bit_size = Int64(skipsize/(3*N))

    if Int(bit_size) == 4
        dtype = Float32
        #println("Reading single precision snapshot")
    elseif Int(bit_size) == 8
        dtype = Float64
        println("Reading double precision snapshot")
    else
        println("read error! neither 32 nor 64 bits data!")
        return -1
    end

    seek(f, 288)

    # set up dictionaries for particles
    for i in 1:length(data["Header"]["PartTypes"])
        if data["Header"]["npart"][i] != 0
            data[data["Header"]["PartTypes"][i]] = Dict()
        end
        if data["Header"]["massarr"][i] != 0
           data[data["Header"]["PartTypes"][i]]["MASS"] = dtype.(data["Header"]["massarr"][i] * ones(data["Header"]["npart"][i]))
        end
    end

    blockname = "POS"


    while (blockname != "INFO") #& (eof(f) != true)

        #println(blockname)
        #println(typeof(blockname))
        # skip identifiers
        p = position(f)

        seek(f, p+8)

        if eof(f) == true
            break
        end

        p, data = read_block(p, data, dtype, blockname, f, bit_size)


        seek(f,p+4)

        if eof(f) == true
            break
        end

        p = position(f)
        seek(f, p+4)
        # read blockname

        name = Char.(read!(f, Array{Int8,1}(undef,4)))
        blockname = String(name)

        blockname = String(strip(blockname))

        p = position(f)
        #seek(f, p+8)

    end


    close(f)

    return data


end

function read_block_with_info(f::IOStream, d::Dict{Any,Any}, info::Info_Line)

    parttypes = ["PartType0", "PartType1", "PartType2",
                 "PartType3", "PartType4", "PartType5"]

    for i ∈ 1:length(info.is_present)

        if info.is_present[i] == Int32(1)
            d[parttypes[i]][info.block_name] = copy(transpose(
                                         read!(f, Array{info.data_type,2}(undef,(info.n_dim,d["Header"]["npart"][i])))))
        end

    end

    return d

end

function read_block(p::Int64, data::Dict{Any,Any}, dtype::DataType, blockname::String,
                    f::IOStream, bit_size::Int64)



    file_curr = @__FILE__
    path_curr = file_curr[1:end-17]

    include(path_curr * "part_specific_fields.jl")

    N = sum(data["Header"]["npart"])

    skipsize = read(f, Int32)

    if blockname == "MASS"
        for i ∈ 1:length(data["Header"]["PartTypes"])

            n = data["Header"]["npart"][i]

            if data["Header"]["npart"][i] != Int32(0)

                if data["Header"]["massarr"][i] == Int32(0)
                    data[data["Header"]["PartTypes"][i]][blockname] = copy(transpose(read!(f, Array{dtype,2}(undef,(1,n)))))
                else
                    data[data["Header"]["PartTypes"][i]][blockname] =  Float32.(data["Header"]["massarr"][i] .* ones(n))
                end

            end
        end
    elseif blockname != "ID"

        gas_block = any(x->x==blockname,gas_arr)

        if gas_block == true

            #println("Reading gas block")

            n = Int64(data["Header"]["npart"][1])

            dim = Int(skipsize/(n*bit_size))

            #println("nr. of gas-particles: ", n)
            #println("reading dimensions: ", skipsize/(n*bit_size))

            data[data["Header"]["PartTypes"][1]][blockname] = copy(transpose(read!(f, Array{dtype,2}(undef,(dim,n)))))

        # check if bh
        elseif (any(x->x==blockname,bh_arr)) == true

            n = Int64(data["Header"]["npart"][6])

            #println("nr. of bh-particles: ", n)

            dim = Int(skipsize/(n*bit_size))

            #println("reading dimensions: ", skipsize/(n*bit_size))

            if blockname != "BHPC"
                data[data["Header"]["PartTypes"][6]][blockname] = copy(transpose(read!(f, Array{dtype,2}(undef,(dim,n)))))
            else
                data[data["Header"]["PartTypes"][6]][blockname] = Int.(copy(transpose(read!(f, Array{UInt32,2}(undef,(1,n))))))
            end

        elseif blockname == gas_bh

            n = Int64(data["Header"]["npart"][1])
            data[data["Header"]["PartTypes"][1]][blockname] = Int.(copy(transpose(read!(f, Array{UInt32,2}(undef,(1,n))))))

            n = Int64(data["Header"]["npart"][6])
            data[data["Header"]["PartTypes"][6]][blockname] = Int.(copy(transpose(read!(f, Array{UInt32,2}(undef,(1,n))))))

        else
            # read Blocks
            #println("Blockname: ", blockname)

            for i in 1:length(data["Header"]["PartTypes"])

                dim = (skipsize/(N*bit_size))

                dim = Int64(trunc(dim))

                if data["Header"]["npart"][i] != 0

                    n = Int64(data["Header"]["npart"][i])

                    data[data["Header"]["PartTypes"][i]][blockname] = copy(transpose(read!(f, Array{dtype,2}(undef,(dim,n)))))

                end

            end
        end

    else

        for i in 1:length(data["Header"]["PartTypes"])

            if data["Header"]["npart"][i] != 0

                n = Int64(data["Header"]["npart"][i])

                data[data["Header"]["PartTypes"][i]][blockname] = copy(transpose(read!(f, Array{UInt32,2}(undef,(1,n)))))

            end

        end

    end

    p = position(f)

    return p, data
end



function snap_1_d(filename::String, data::Dict{Any,Any})

    println("wrong function call")
    f = open(filename)

    seek(f, 264)

    N = sum(data["Header"]["npart"])
    skipsize = read(f, Int32)
    bit_size = Int64(skipsize/(3*N))

    if Int(bit_size) == 4
        dtype = Float32
        println("Reading single precision snapshot")
    elseif Int(bit_size) == 8
        dtype = Float64
        println("Reading double precision snapshot")
    else
        println("read error! neither 32 nor 64 bits data!")
        return -1
    end

    seek(f, 288)

    # set up dictionaries for particles
    for i in 1:length(data["Header"]["PartTypes"])
        if data["Header"]["npart"][i] != 0
            data[data["Header"]["PartTypes"][i]] = Dict()
        end
    end

    blockname = "POS"

    while eof(f) != true

        #println(blockname)
        # skip identifiers
        p = position(f)

        p, data = read_block(p, data, dtype, blockname, f, bit_size)

        seek(f,p+4)

        if eof(f) == true
            break
        end

        # read blockname
        letters = read!(f, read!(f, Array{Int8,1}(undef,4)))
        block = Char.(letters)
        block = string.(block)

        blockname = lowercase(strip(block[1] * block[2] * block[3] * block[4]))
    end


    close(f)

end

function snap_1_d_old(filename, data)

    println("correct function call")
    f = open(filename)

    seek(f, 264)

    N = sum(data["Header"]["npart"])
    skipsize = read(f, Int32)
    bit_size = skipsize/(3*N)

    if Int(bit_size) == 4

    elseif Int(bit_size) == 8

    else
        println("read error! neither 32 nor 64 bits data!")
        return -1
    end

    # read positions
    for i in 1:length(data["Header"]["PartTypes"])

        if data["Header"]["npart"][i] != 0

            data[data["Header"]["PartTypes"][i]] = Dict()
            N = Int64(data["Header"]["npart"][i])
            #dummy = read(f, Float32, (3,N))

            if Int(bit_size) == 4
                data[data["Header"]["PartTypes"][i]]["POS"] = copy(transpose(read!(f, Array{Float32,2}(undef,(3,N)))))
            else
                data[data["Header"]["PartTypes"][i]]["POS"] = copy(transpose(read!(f, Array{Float64,2}(undef,(3,N)))))
            end


        end

    end

    p = position(f)

    # skip identifiers
    seek(f, p+8)

    # read Velocities
    for i in 1:length(data["Header"]["PartTypes"])

        if data["Header"]["npart"][i] != 0

            N = Int64(data["Header"]["npart"][i])

            if Int(bit_size) == 4
                data[data["Header"]["PartTypes"][i]]["VEL"] = copy(transpose(read!(f, Array{Float32,2}(undef,(3,N)))))
            else
                data[data["Header"]["PartTypes"][i]]["VEL"] = copy(transpose(read!(f, Array{Float64,2}(undef,(3,N)))))
            end

        end

    end

    #dummy = 0
    #gc()

    p = position(f)

    # skip identifiers
    seek(f, p+8)

    # Read IDs
    for i in 1:length(data["Header"]["PartTypes"])

        if data["Header"]["npart"][i] != 0

            N = Int64(data["Header"]["npart"][i])
            data[data["Header"]["PartTypes"][i]]["ID"] = copy(transpose(read!(f, Array{UInt32,1}(undef,N))))

        end

    end

    p = position(f)

    # skip identifiers
    seek(f, p+8)

    # Read IDs
    for i in 1:length(data["Header"]["PartTypes"])

        if data["Header"]["npart"][i] != 0

            if data["Header"]["massarr"][i] == 0

                N = Int64(data["Header"]["npart"][i])

                if Int(bit_size) == 4
                    data[data["Header"]["PartTypes"][i]]["MASS"] = copy(transpose(read!(f, Array{Float32,1}(undef,N))))
                else
                    data[data["Header"]["PartTypes"][i]]["MASS"] = copy(transpose(read!(f, Array{Float64,1}(undef,N))))
                end

            else

                N = Int64(data["Header"]["npart"][i])
                data[data["Header"]["PartTypes"][i]]["MASS"] = zeros(N)
                data[data["Header"]["PartTypes"][i]]["MASS"] .= data["Header"]["massarr"][i]

            end

        end

    end

    # if data["Header"]["npart"][1] != 0
    #
    #     # skip identifiers
    #     p = position(f)
    #     seek(f, p+8)
    #
    #     N = Int64(data["Header"]["npart"][0])
    #
    #     # read U
    #     if Int(bit_size) == 4
    #         data["PartType0"]["InternalEnergy"] = read(f, Float32, N)
    #     else
    #         data["PartType0"]["InternalEnergy"] = read(f, Float64, N)
    #     end
    #
    #
    #     # skip identifiers
    #     p = position(f)
    #     seek(f, p+8)
    #
    #     # Read Density
    #     if Int(bit_size) == 4
    #         data["PartType0"]["Density"] = read(f, Float32, N)
    #     else
    #         data["PartType0"]["Density"] = read(f, Float64, N)
    #     end
    #
    #     # skip identifiers
    #     p = position(f)
    #     seek(f, p+8)
    #
    #     # Read SmoothingLength
    #     if Int(bit_size) == 4
    #         data["PartType0"]["SmoothingLength"] = read(f, Float32, N)
    #     else
    #         data["PartType0"]["SmoothingLength"] = read(f, Float64, N)
    #     end
    #
    # end

    close(f)

    return data


end

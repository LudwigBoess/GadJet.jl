
"""
    read_block_by_name(filename::String, blockname::String;
                                info::Info_Line=Info_Line(),
                                parttype::Int64=-1)

Reads a block in a snapshot with given name. Names are case sensitive.

# Examples
```jldoctest
julia> pos_info = Info_Line("POS", Float32, 1, [1, 1, 1, 1, 1, 1])
[...]
julia> gas_pos = read_block_by_name(filename, "POS", info=pos_info, parttype=0)
[...]
```
"""
function read_block_by_name(filename::String, blockname::String;
                            info::Info_Line=Info_Line(),
                            parttype::Int64=-1)


    # read header - super fast and needed for flexibility
    h = head_to_obj(filename)

    if parttype == -1
        parttypes = ["PartType0", "PartType1", "PartType2",
                     "PartType3", "PartType4", "PartType5"]
        d = Dict()
    end # if parttype

    blockname = strip(blockname)

    if info.block_name == ""
        info = read_info(filename)
        if info == 1
            if blockname == "MASS"
                info = Info_Line("MASS", Float32, 1, [0, 0, 0, 0, 0, 0])
            else
                error("No Info block in snapshot! Supply Info_Line type!")
            end
        else # info != 1
            for i ∈ 1:length(info)
                if info[i].block_name == blockname
                    info = info[i]
                    break
                end # if block found
            end # loop over info
            if isa(info, Array)
                if (blockname == "MASS")
                    info = Info_Line("MASS", Float32, 1, [0, 0, 0, 0, 0, 0])
                else
                    error("Block not present!")
                end
            end
        end # info == 1
    end

    f = open(filename)

    first_block = read(f, Int32)

    if first_block != Int32(8)
        error("Only snapshots of format 2 can be read by name!")
        return
    end

    while eof(f) != true
        name = Char.(read!(f, Array{Int8,1}(undef,4)))
        blname = String(name)
        blname = strip(blname)

        if blname == blockname
            break
        end
        p = position(f)
        seek(f,p+8)
        skipsize = read(f, Int32)
        p = position(f)
        seek(f,p+skipsize+8)
    end

    if eof(f) != true
        p = position(f)
        seek(f,p+8)

        blocksize = read(f,Int32)

        for i ∈ 1:length(info.is_present)
            p = position(f)

            if info.is_present[i] == Int32(1)
                if parttype == -1
                    d[parttypes[i]] = Dict()
                    d[parttypes[i]][blockname] = copy(transpose(
                                                        read!(f, Array{info.data_type,2}(undef,(info.n_dim,h.npart[i])))))
                else
                    if i == (parttype+1)

                        block = copy(transpose(
                                        read!(f, Array{info.data_type,2}(undef,(info.n_dim,h.npart[i])))))
                        close(f)
                        return block
                    else
                        seek(f, p + ( sizeof(info.data_type)*info.n_dim*h.npart[i] ))
                    end # if i == (parttype+1)
                end # parttype == -1
            end # info.is_present[i] == Int32(1)
        end # i ∈ 1:length(info.is_present)

        # fill in mass from header
        if blockname == "MASS"
            for i ∈ 1:length(info.is_present)
                if info.is_present[i] == Int32(0)
                    if parttype == -1
                        d[parttypes[i]] = Dict()
                        d[parttypes[i]][blockname] = Array{info.data_type,2}(undef,(info.n_dim,h.npart[i]))
                        d[parttypes[i]][blockname] .= h.massarr[i]
                    end # parttype == -1
                end # info block
            end # loop over i
        end
        close(f)

        return d

    else # blockname not found
        if blockname != "MASS"
            error("Block not present!")
        else
            if parttype == -1
                for i = 1:length(parttypes)
                    d[parttypes[i]] = Dict()
                    d[parttypes[i]][blockname] = Array{info.data_type,2}(undef,(info.n_dim,h.npart[i]))
                    d[parttypes[i]][blockname] .= h.massarr[i]
                end
            else
                block = Array{info.data_type,2}(undef,(info.n_dim,h.npart[parttype+1]))
                block .= h.massarr[parttype+1]
                return block
            end # parttype == -1
        end
    end # eof(f) != true

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


function read_block_with_offset(filename::String, data_old, pos0::Integer, info::Info_Line,
                                offset::Integer, offset_key, n_to_read::Integer, part_per_key )

    # open the file
    f = open(filename)

    # number of bits in data_type
    len = sizeof(info.data_type) * info.n_dim

    # jump to position of particle type in relevant block
    seek(f, pos0+offset*len)

    # store position in file
    p = position(f)

    # allocate array to store data
    data =  Array{info.data_type,2}(undef, (n_to_read, info.n_dim))

    n_read = 1
    n_this_key = 0

    for i = 1:length(offset_key)

        # jump to start of key
        seek(f, p + len*offset_key[i])
        n_this_key += part_per_key[i]

        data[n_read:n_this_key, :] = copy(transpose(read!(f,
                    Array{info.data_type,2}(undef,(info.n_dim,part_per_key[i])))))

        n_read += part_per_key[i]

    end # for

    # close the file
    close(f)

    return [data_old; data]
end

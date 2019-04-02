#__precompile__()

module GadJet


    include(joinpath(dirname(@__FILE__), "read_snapshot", "gadget_types.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "dict_functions.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "obj_functions.jl"))
    include(joinpath(dirname(@__FILE__), "write_snapshot", "write_snap.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_to_grid.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_types.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "kernels.jl"))
    include(joinpath(dirname(@__FILE__), "unit_conversion", "unit_types.jl"))


    export Header, Info_Line,       # types
           head_to_dict,
           snap_to_dict,
           head_to_obj,
           print_blocks,
           read_info,
           read_block_by_name,
           readnew,                 # similar to readnew.pro by Klaus Dolag
           write_header,
           write_block,
           sphAdaptiveMapping,
           sphCenterMapping,
           Cubic,
           Quintic,
           WendlandC4,
           WendlandC6,
           mappingParameters,
           GadgetUnitFactors



    # returns an array of present blocknames
    function print_blocks(filename)

        f = open(filename)
        blocksize = read(f, Int32)

        if blocksize != 8
            return "Not possible - use snap_format 2!"
        end

        p = position(f)
        seek(f,4)

        blocks = []

        while eof(f) != true

            name = Char.(read!(f, Array{Int8,1}(undef,4)))
            blockname = String(name)

            blockname = strip(blockname)

            push!(blocks, blockname)

            p = position(f)
            seek(f,p+8)

            skipsize = read(f, Int32)

            p = position(f)
            seek(f,p+skipsize+8)

        end

        println("Found blocks: ")
        for block ∈ blocks
            println(block)
        end

        return blocks

    end

    function read_info(filename; verbose::Bool=false)

        f = open(filename)
        seek(f,4)

        while eof(f) != true

            name = Char.(read!(f, Array{Int8,1}(undef,4)))
            blockname = String(name)

            p = position(f)
            seek(f,p+8)

            skipsize = read(f, Int32)

            if blockname == "INFO"
                n_blocks = Int(skipsize/40)
                arr_info = Array{Info_Line,1}(undef,n_blocks)

                for i = 1:n_blocks
                    arr_info[i] = read_info_line(f)
                end # for

                close(f)

                if verbose == true
                    println("Found Info block.\nEntries are:\n\n")
                    for i ∈ 1:length(arr_info)
                        println(i, " - ", arr_info[i].block_name)
                    end # for
                end # verbose

                return arr_info

            else
                p = position(f)
                seek(f,p+skipsize+8)
            end # if blockname == "INFO"

        end # while eof(f) != true

        println("No info block present!")

        return 1

    end

    function read_info_line(f)

        name = Char.(read!(f, Array{Int8,1}(undef,4)))
        blockname = String(name)

        block_name = strip(blockname)

        letters = read!(f, Array{Int8,1}(undef,8))
        letters = Char.(letters)
        letters = string.(letters)

        data_type = ""
        for i = 1:length(letters)
            data_type *= letters[i]
        end

        data_type = lowercase(strip(data_type))

        if data_type == "float" ||
           data_type == "floatn"
                dt = Float32
        elseif data_type == "long"
                dt = UInt32
        elseif data_type == "double" ||
               data_type == "doublen"
                dt = Float64
        end

        n_dim = read(f,Int32)

        is_present = read!(f, Array{Int32,1}(undef,6))

        in_l = Info_Line(block_name, dt, n_dim, is_present)

        return in_l
    end


    function read_block_by_name(filename::String, blockname::String;
                                info::Info_Line=Info_Line(),
                                npart::Vector{Int32}=Int32.(zeros(6)),
                                parttype::Int64=-1)

        
        if parttype == -1
            parttypes = ["PartType0", "PartType1", "PartType2",
                         "PartType3", "PartType4", "PartType5"]
        end # if parttype

        blockname = uppercase(strip(blockname))

        if info.block_name == ""
            info = read_info(filename)
            if info == 1
                println("Error! No Info block in snapshot! Supply Info_Line type!")
                return 1
            end
            for i ∈ 1:length(info)
                if info[i].block_name == blockname
                    info = info[i]
                    break
                end
            end
        end

        if npart == Int32.(zeros(6))
            h = head_to_obj(filename)
            npart = h.npart
        end

        f = open(filename)

        first_block = read(f, Int32)

        if first_block != Int32(8)
            return "Only snapshots of format 2 can be read by name!"
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

            d = Dict()

            for i ∈ 1:length(info.is_present)
                p = position(f)

                if info.is_present[i] == Int32(1)
                    if parttype == -1
                        d[parttypes[i]] = Dict()
                        d[parttypes[i]][blockname] = copy(transpose(
                                                            read!(f, Array{info.data_type,2}(undef,(info.n_dim,npart[i])))))
                    else
                        if i == (parttype+1)

                            block = copy(transpose(
                                            read!(f, Array{info.data_type,2}(undef,(info.n_dim,npart[i])))))
                            close(f)
                            return block
                        else
                            seek(f, p + ( sizeof(info.data_type)*info.n_dim*npart[i] ))
                        end # if i == (parttype+1)
                    end # parttype == -1
                end # info.is_present[i] == Int32(1)
            end # i ∈ 1:length(info.is_present)


            close(f)

            return d
        else

            return "Block not present!"
        end # eof(f) != true

    end


    function readnew(filename::String, info::Info_Line, npart::Vector{Int32};
                     parttype::Int64=0)

        f = open(filename)

        blockname = strip(info.block_name)

        first_block = read(f, Int32)

        if first_block != Int32(8)
            println("Error! Only snapshots of format 2 can be read by name!")
            return 1
        end

        while eof(f) != true

            name = strip(String(Char.(read!(f, Array{Int8,1}(undef,4)))))

            if name == blockname
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

            d = Dict()

            for i ∈ 1:length(info.is_present)
                p = position(f)
                if (info.is_present[i] == Int32(1)) && (i == Int(parttype + 1))
                    d = copy(transpose(
                            read!(f, Array{info.data_type,2}(undef,(info.n_dim,npart[i])))))
                    close(f)
                    return d
                else
                    seek(f, p + ( sizeof(info.data_type)*info.n_dim*npart[i] ) )
                end

            end


            close(f)
            println("Error! Block not present for particle type!")
            return 1
        else
            println("Block not present!")
            return 1
        end

    end


end

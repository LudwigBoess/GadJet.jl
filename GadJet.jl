#__precompile__()

module GadJet


    include(joinpath(dirname(@__FILE__), "read_snapshot", "gadget_types.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "dict_functions.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "obj_functions.jl"))
    include(joinpath(dirname(@__FILE__), "write_snapshot", "write_snap.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_to_grid.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_types.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "kernels.jl"))


    export Header, Info_Line,       # types
           head_to_dict,
           snap_to_dict,
           head_to_obj,
           print_blocks,
           get_info,
           read_block_by_name,
           write_header,
           write_block,
           sphAdaptiveMapping,
          # sph_to_grid,
           Cubic,
           Quintic,
           WendlandC4,
           WendlandC6,
           mappingParameters



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

    function get_info(filename)

        f = open(filename)

        seekend(f)
        p = position(seekend(f))

        seek(f, p-4)
        blocksize = read(f, Int32)

        # one info line is 40 bytes long
        if blocksize%40 == 0

            seek(f, p-blocksize-4)
            #blocksize = read(f,Int32)
            n_blocks = Int(blocksize/40)
            arr_info = Array{Info_Line,1}(undef,n_blocks)

            for i = 1:n_blocks
                arr_info[i] = get_info_line(f)
            end

            blocksize = read(f, Int32)
            close(f)

            return arr_info
        else

            return "No info block present!"
        end

    end

    function get_info_line(f)

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

    function read_block_by_name(filename::String, blockname::String,
                                info::Info_Line, npart::Vector{Int32})

        parttypes = ["PartType0", "PartType1", "PartType2",
                     "PartType3", "PartType4", "PartType5"]

        f = open(filename)

        blockname = strip(blockname)

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

                if info.is_present[i] == Int32(1)
                    d[parttypes[i]] = Dict()
                    d[parttypes[i]][blockname] = copy(transpose(
                                                        read!(f, Array{info.data_type,2}(undef,(info.n_dim,npart[i])))))
                end

            end


            close(f)

            return d
        else

            return "Block not present!"
        end

    end


end

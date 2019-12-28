#__precompile__()

module GadJet


    include(joinpath(dirname(@__FILE__), "read_snapshot", "gadget_types.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "dict_functions.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "obj_functions.jl"))
    include(joinpath(dirname(@__FILE__), "write_snapshot", "write_snap.jl"))
    # sph to grid mapping internal module
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "kernels.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_to_grid.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_types.jl"))
    # sph to grid mapping with Smac
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "smac1_utility.jl"))
    # sph to grid mapping with Smac
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "smac2_utility.jl"))
    # unit conversion
    include(joinpath(dirname(@__FILE__), "unit_conversion", "unit_types.jl"))
    # old riemann solver
    #include(joinpath(dirname(@__FILE__), "riemann_solvers", "riemann_solver.jl"))

    # test for new riemann solvers
    include(joinpath(dirname(@__FILE__), "riemann_solvers", "setup_riemann_parameters.jl"))
    #include(joinpath(dirname(@__FILE__), "riemann_solvers", "cr_sod_shock.jl"))

    include(joinpath(dirname(@__FILE__), "bp_cr_utility", "cr_datatypes.jl"))
    include(joinpath(dirname(@__FILE__), "bp_cr_utility", "analysis_functions.jl"))


    export Header, Info_Line,       # types
           head_to_dict,
           snap_to_dict,
           head_to_obj,
           print_blocks,
           read_info,
           read_block_by_name,      # similar to readnew.pro by Klaus Dolag
           write_header,
           write_block,

           # Kernels
           Cubic,
           Quintic,
           WendlandC4,
           WendlandC6,
           # internal sph mapping
           mappingParameters,
           sphAdaptiveMapping,
           sphCenterMapping,
           # helper functions and datatypes for Smac
           Smac1ImageInfo,
           read_smac1_binary_image,
           read_smac1_binary_info,
           # helper function for P-Smac2
           write_smac2_par,

           GadgetUnitFactors,
           # old riemann solver
           #RiemannParameters,   # datatype for riemann parameters
           #RiemannSolution,     # datatype for riemann solution
           #solveHydroShock,      # function that solves a standard sod shock

           # Test for different riemann solvers
           RiemannParameters,    # helper function to set up solution
           solve,                # overloaded function to solve riemann problems

           # datatypes and helper functions for BP_REAL_CRs
           CRShockData,          # datatype to analyse single shocked particle
           readSingleCRShockDataFromOutputFile, # as the name says
           CRMomentumDistributionConfig, # config parameters for momentum distribution function
           CRMomentumDistribution,
           getCRMomentumDistributionFromPartID, # function to get distribution function
           calculateCREnergyInCGS,
           calculateCRNumber




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
                    println("Found Info block.\nEntries are:\n")
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
                dt = Int32
        elseif data_type == "llong"
                dt = Int64
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

        blockname = strip(blockname)

        if info.block_name == ""
            info = read_info(filename)
            if info == 1
                error("No Info block in snapshot! Supply Info_Line type!")
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

end
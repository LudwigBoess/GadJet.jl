# include("cr_datatypes.jl")
# file_curr = @__FILE__
# path_gadjet = file_curr[1:end-35]
# include(path_gadjet * "GadJet.jl")

function readSingleCRShockDataFromOutputFile(file)

        # read file into memory
        f = open(file)
        lines = readlines(f)
        close(f)

        # filter only relevant lines. Output of every line starts with "CR DATA".
        g = occursin.("CR DATA", lines)
        lines = lines[g]

        # get number of output lines
        N = length(lines)

        # init datatype for analyze
        cr = CRShockData(N)

        for i ∈ 1:N
                line_content            = split(lines[i])
                cr.dt[i]                = parse(Float64,line_content[3])
                cr.Mach[i]              = parse(Float64,line_content[4])
                cr.Shock_Speed[i]       = parse(Float64,line_content[5])
                cr.Shock_Compress[i]    = parse(Float64,line_content[6])
                cr.Shock_Energy_In[i]   = parse(Float64,line_content[7])
                cr.Shock_Energy_Real[i] = parse(Float64,line_content[8])
                cr.Energy_P[i]          = parse(Float64,line_content[9])
                cr.Energy_e[i]          = parse(Float64,line_content[10])
        end

        return cr
end

function getCRMomentumDistributionFromPartID(snap_file::String, ID::Int64;
                                             pmin::Float64=10.0, pmax::Float64=1.0e7,
                                             Nbins::Int64=0)


    ymin = 1.e-20
    h = head_to_obj(snap_file)
    info = read_info(snap_file)

    if info == 1
        if Nbins == 0
            println("Can't read spectrum! No info block present!\nSupply number of momentum bins to proceed!")
            return 1
        else
            info = Array{Info_Line,1}(undef,7)
            info[1] = Info_Line("ID", UInt32, Int32(1), [1, 0, 0, 0, 0, 0])
            info[2] = Info_Line("CRpN", Float32, Int32(24), [1, 0, 0, 0, 0, 0])
            info[3] = Info_Line("CRpS", Float32, Int32(24), [1, 0, 0, 0, 0, 0])
            info[4] = Info_Line("CRpC", Float32, Int32(1), [1, 0, 0, 0, 0, 0])
            info[5] = Info_Line("CReN", Float32, Int32(24), [1, 0, 0, 0, 0, 0])
            info[6] = Info_Line("CReS", Float32, Int32(24), [1, 0, 0, 0, 0, 0])
            info[7] = Info_Line("CReC", Float32, Int32(1), [1, 0, 0, 0, 0, 0])
        end
    end

    id = read_block_by_name(snap_file, "ID",
                            info=info[getfield.(info, :block_name) .== "ID"][1],
                            npart=h.npart, parttype=0)

    part = findfirst( id .== UInt32(ID) )

    # protons
    CRpN = read_block_by_name(snap_file, "CRpN",
                            info=info[getfield.(info, :block_name) .== "CRpN"][1],
                            npart=h.npart, parttype=0)[part,:]
    CRpS = read_block_by_name(snap_file, "CRpS",
                            info=info[getfield.(info, :block_name) .== "CRpS"][1],
                            npart=h.npart, parttype=0)[part,:]
    CRpC = read_block_by_name(snap_file, "CRpC",
                            info=info[getfield.(info, :block_name) .== "CRpC"][1],
                            npart=h.npart, parttype=0)[part]

    # electrons
    CReN = read_block_by_name(snap_file, "CReN",
                            info=info[getfield.(info, :block_name) .== "CReN"][1],
                            npart=h.npart, parttype=0)[part,:]
    CReS = read_block_by_name(snap_file, "CReS",
                            info=info[getfield.(info, :block_name) .== "CReS"][1],
                            npart=h.npart, parttype=0)[part,:]
    CReC = read_block_by_name(snap_file, "CReC",
                            info=info[getfield.(info, :block_name) .== "CReC"][1],
                            npart=h.npart, parttype=0)[part]

    Nbins = length(CRpS)
    par = CRMomentumDistributionConfig(pmin, pmax, Nbins)

    cr = CRMomentumDistribution(par.Nbins)

    # get zeroth bin
    cr.CRp_bound[1] = pmin
    cr.CRp_dis[1] = CRpN[1]

    cr.CRe_bound[1] = pmin
    cr.CRe_dis[1] = CReN[1]

    # all other bins
    j = 2
    for i = 1:Nbins-1

        # upper boundary of bin
        cr.CRp_bound[j] = pmin * 10.0^(par.bin_width*i)
        cr.CRp_dis[j] = CRpN[i] * ( cr.CRp_bound[j]/cr.CRp_bound[j-1])^(-CRpS[i])
        if cr.CRp_bound[j] > CRpC/par.mc_p
            cr.CRp_bound[j] = CRpC/par.mc_p
        end

        cr.CRe_bound[j] = pmin * 10.0^(par.bin_width*i)
        cr.CRe_dis[j] = CReN[i] * ( cr.CRe_bound[j]/cr.CRe_bound[j-1])^(-CReS[i])
        if cr.CRe_bound[j] > CReC/par.mc_e
            cr.CRe_bound[j] = CReC/par.mc_e
        end

        # lower bound of next bin
        cr.CRp_bound[j+1] = cr.CRp_bound[j]
        cr.CRp_dis[j+1] = CRpN[i+1]
        if cr.CRp_bound[j] == CRpC/par.mc_p
            cr.CRp_bound[j+1] = CRpC/par.mc_p
        end

        cr.CRe_bound[j+1] = cr.CRe_bound[j]
        cr.CRe_dis[j+1] = CReN[i+1]
        if cr.CRe_bound[j] == CReC/par.mc_e
            cr.CRe_bound[j+1] = CReC/par.mc_e
        end

        j += 2

    end

    # last boundary
    cr.CRp_bound[j] = pmax
    if cr.CRp_bound[j-1] < CRpC/par.mc_p
        cr.CRp_bound[j] = CRpC/par.mc_p
    end
    cr.CRp_dis[j] = CRpN[Nbins] * ( cr.CRp_bound[j]/cr.CRp_bound[j-1])^(-CRpS[Nbins])
    cr.CRp_bound[j+1] = cr.CRp_bound[j]
    cr.CRp_dis[j] = ymin

    cr.CRe_bound[j] = pmax
    if cr.CRe_bound[j-1] < CReC/par.mc_e
        cr.CRe_bound[j] = CReC/par.mc_e
    end
    cr.CRe_dis[j] = CReN[Nbins] * ( cr.CRe_bound[j]/cr.CRe_bound[j-1])^(-CReS[Nbins])
    cr.CRe_bound[j+1] = cr.CRe_bound[j]
    cr.CRe_dis[j] = ymin * 1.e-2

    return cr
end

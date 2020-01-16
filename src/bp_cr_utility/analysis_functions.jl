
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
                            parttype=0)

    part = findfirst( id .== UInt32(ID) )[1]

    # protons
    CRpN = read_block_by_name(snap_file, "CRpN",
                            info=info[getfield.(info, :block_name) .== "CRpN"][1],
                            parttype=0)[part,:]
    CRpS = read_block_by_name(snap_file, "CRpS",
                            info=info[getfield.(info, :block_name) .== "CRpS"][1],
                            parttype=0)[part,:]
    CRpC = read_block_by_name(snap_file, "CRpC",
                            info=info[getfield.(info, :block_name) .== "CRpC"][1],
                            parttype=0)[part]

    # electrons
    CReN = read_block_by_name(snap_file, "CReN",
                            info=info[getfield.(info, :block_name) .== "CReN"][1],
                            parttype=0)[part,:]
    CReS = read_block_by_name(snap_file, "CReS",
                            info=info[getfield.(info, :block_name) .== "CReS"][1],
                            parttype=0)[part,:]
    CReC = read_block_by_name(snap_file, "CReC",
                            info=info[getfield.(info, :block_name) .== "CReC"][1],
                            parttype=0)[part]

    Nbins = length(CRpS)

    par = CRMomentumDistributionConfig(pmin, pmax, Nbins)

    cr = CRMomentumDistribution(par.Nbins)

    # compensate for io - needs to be converted to Float64!
    CRpN = 1.e20 .* Float64.(CRpN)
    CReN = 1.e20 .* Float64.(CReN)

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
    #cr.CRp_dis[j] = ymin

    cr.CRe_bound[j] = pmax
    if cr.CRe_bound[j-1] < CReC/par.mc_e
        cr.CRe_bound[j] = CReC/par.mc_e
    end
    cr.CRe_dis[j] = CReN[Nbins] * ( cr.CRe_bound[j]/cr.CRe_bound[j-1])^(-CReS[Nbins])
    cr.CRe_bound[j+1] = cr.CRe_bound[j]
    #cr.CRe_dis[j] = ymin * 1.e-2

    return cr
end


function energy_integral(bound_low, bound_up, norm, slope, ρ)
    # energy integral (eq.21 M01)

    cnst_c = 2.9979e10
    slope_soft = 1.e-6

    en = 4.0 * π * cnst_c * norm * bound_low^4 /  ρ
    energy = en * ( (bound_up/bound_low)^(4.0 - slope ) - 1.0 ) / ( 4.0 - slope )

    if ( 4.0 - slope_soft ) < slope < ( 4.0 + slope_soft )
        slope_var = (slope - 4.0)/slope_soft
        energy2 = en * log10(bound_up/bound_low)
        if slope_var != 0.0
            energy = energy * slope_var + energy2 * ( 1.0 - slope_var )
        else
            energy = energy2
        end
    end

    return energy
end


function calculateCREnergyInCGS(CR_N, CR_S, CR_Cut, ρ; pmin=10.0, pmax=1.e7, mc=0.0, SelectBin=-1, units::GadgetPhysicalUnits=GadgetPhysicalUnits(), verbose=true)
    # calculates engergy per bin in cgs units.

    if mc == 0.0
        println("Error! No mc specified! Please provide m_particle * c_light to continue!")
        return 1
    end

    CR_E = 0.0
    cnst_c = 2.9979e10
    ρ *= units.m_unit / units.l_unit^3  # convert rho in cgs units

    CR_N .*= 1.e20

    Nbins = length(CR_N)
    bin_width = log10(pmax/pmin)/Nbins
    bin = 10.0.^collect(log10(pmin):bin_width:log10(pmax)) .* mc
    above_cut = findall( bin .> CR_Cut/mc )
    if length(above_cut) > 0
        bin[above_cut] .= CR_Cut/mc
    end

    if SelectBin == -1
        CR_E = zeros(length(CR_N))
        for i = 1:length(CR_N)
            CR_E[i] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
        end
    elseif length(SelectBin) == 1
        CR_E = energy_integral(bin[SelectBin], bin[SelectBin+1], CR_N[SelectBin], CR_S[SelectBin], ρ)
    elseif length(SelectBin) > 1
        CR_E = zeros(length(SelectBin))
        if i ∈ SelectBin
            CR_E[i-SelectBin[1]+1] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
        end
    end
    return CR_E
end


function density_integral(bound_low, bound_up, norm, slope, ρ)
    # density integral (eq. 9 M01)

    slope_soft = 1.e-6

    nb = 4.0 * π * norm * bound_low^3 / ρ
    density = nb * ( (bound_up/bound_low)^(3.0 - slope ) - 1.0 ) / ( 3.0 - slope )

    if ( 3.0 - slope_soft ) < slope < ( 3.0 + slope_soft )
        slope_var = (slope - 3.0)/slope_soft
        density2 = nb * log10(bound_up/bound_low)
        if slope_var != 0.0
            density = density * slope_var + density2 * ( 1.0 - slope_var )
        else
            density = density2
        end
    end

    return density
end


function calculateCRNumber(CR_N, CR_S, CR_Cut, ρ; pmin=10.0, pmax=1.e7, mc, SelectBin=-1, units::GadgetPhysicalUnits=GadgetPhysicalUnits(), verbose=true)
    # calculates engergy per bin in cgs units.

    CR_Num = 0.0
    cnst_c = 2.9979e10
    ρ *= units.m_unit / units.l_unit^3  # convert rho in cgs units

    CR_N .*= 1.e20  # compensate for io

    Nbins = length(CR_N)
    bin_width = log10(pmax/pmin)/Nbins
    bins = 10.0.^collect(log10(pmin):bin_width:log10(pmax)) .* mc
    above_cut = findall( bins .> CR_Cut/mc )

    if length(above_cut) > 0
        bins[above_cut] .= CR_Cut/mc
    end

    if SelectBin == -1
        CR_Num = zeros(length(CR_N))
        for i = 1:length(CR_N)
            CR_Num[i] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
        end
    elseif length(SelectBin) == 1
        CR_Num = energy_integral(bin[SelectBin], bin[SelectBin+1], CR_N[SelectBin], CR_S[SelectBin], ρ)
    elseif length(SelectBin) > 1
        CR_Num = zeros(length(SelectBin))
        if i ∈ SelectBin
            CR_Num[i-SelectBin[1]+1] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
        end
    end
    return CR_Num
end

"""
        INFO HERE
"""

"""
                Debug data from shock
"""
mutable struct CRShockData
   """   Data structure to analyze a single shocked particle.
         To be used with compile option: LMB_CR_DEBUG_TRACK_PARTICLE
   """
   dt::Vector{Float64}
   Mach::Vector{Float64}
   Shock_Speed::Vector{Float64}
   Shock_Compress::Vector{Float64}
   Shock_Energy_In::Vector{Float64}
   Shock_Energy_Real::Vector{Float64}
   Energy_P::Vector{Float64}
   Energy_e::Vector{Float64}

   function CRShockData(N::Int64=0)

       new(zeros(N), zeros(N), zeros(N), zeros(N),
           zeros(N), zeros(N), zeros(N), zeros(N))

   end
end




"""
                Distribution spectrum in momentum space
"""

struct CRMomentumDistributionConfig
   """   Config object to obtain distribution spectrum in momentum space.
   """
   pmin::Float64
   pmax::Float64
   Nbins::Int64
   bin_width::Float64
   mp::Float64
   me::Float64
   c::Float64
   mc_e::Float64
   mc_p::Float64

   function CRMomentumDistributionConfig(pmin::Float64=0.0, pmax::Float64=0.0, Nbins::Int64=24)

       CNST_ME = 9.1095e-28
       CNST_MP = 1.6726e-24
       CNST_C = 2.9979e10
       MCe = CNST_ME * CNST_C
       MCp = CNST_MP * CNST_C
       bin_width = log10(pmax/pmin)/Nbins

       new(pmin, pmax, Nbins, bin_width, CNST_MP, CNST_ME, CNST_C, MCe, MCp)

   end
end

mutable struct CRMomentumDistribution

    CRp_dis::Vector{Float64}
    CRp_bound::Vector{Float64}
    CRe_dis::Vector{Float64}
    CRe_bound::Vector{Float64}

    function CRMomentumDistribution(Nbins::Int64)

        CRp_dis = zeros(Int(2*Nbins))
        CRe_dis = zeros(Int(2*Nbins))
        CRp_bound = zeros(Int(2*Nbins + 1))
        CRe_bound = zeros(Int(2*Nbins + 1))

        new( CRp_dis, CRp_bound,
             CRe_dis, CRe_bound )
    end
end

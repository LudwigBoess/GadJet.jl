#using Roots

"""
    Datatypes for IC Parameters and Solution
"""

mutable struct SodCRSolution

    x::Array{Float64,1}
    rho::Array{Float64,1}
    rho3::Float64
    rho4::Float64
    P_tot::Array{Float64,1}
    P_th::Array{Float64,1}
    P_cr::Array{Float64,1}
    P_cr_p::Array{Float64,1}
    P_cr_e::Array{Float64,1}
    P34_tot::Float64
    P3_th::Float64
    P3_cr::Float64
    P4_th::Float64
    P4_cr::Float64
    U_tot::Array{Float64,1}
    U_th::Array{Float64,1}
    E_cr::Array{Float64,1}
    E_cr_p::Array{Float64,1}
    E_cr_e::Array{Float64,1}
    v::Array{Float64,1}
    v34::Float64
    vt::Float64
    vs::Float64
    xs::Float64
    xr::Float64
    Mach::Float64

    function SodCRSolution(x::Array{Float64,1})
        N = length(x)
        new(x,
            zeros(N),   # rho
            0.0,        # rho3
            0.0,        # rho4
            zeros(N),   # P_tot
            zeros(N),   # P_th
            zeros(N),   # P_cr
            zeros(N),   # P_cr_p
            zeros(N),   # P_cr_e
            0.0,        # P34_tot
            0.0,        # P3_th
            0.0,        # P3_cr
            0.0,        # P4_th
            0.0,        # P4_cr
            zeros(N),   # U_tot
            zeros(N),   # U_th
            zeros(N),   # E_cr
            zeros(N),   # E_cr_p
            zeros(N),   # E_cr_e
            zeros(N),   # v
            0.0,        # v34
            0.0,        # vt
            0.0,        # vs
            0.0,        # xs
            0.0,        # xr
            0.0)        # Mach
    end
end

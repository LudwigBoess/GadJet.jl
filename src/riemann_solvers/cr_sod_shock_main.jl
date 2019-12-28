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


"""
        DSA Models
"""

"""
    Kang&Ryu 2007, http://arxiv.org/abs/0704.1521v1
"""
@inline function kr_fitting_function(M::Float64, p::Array{Float64,1})

    mm = M - 1.0
    m2 = M * M

    return ( p[1] + p[2]*mm + p[3]*mm*mm + p[4]*mm*mm*mm + p[5]*mm*mm*mm*mm ) / (m2*m2)
end

function KR07_acc(M::Float64)
    if M <= 2.0
        return 1.96e-3*(M*M - 1.)             # eq. A3
    else
        return kr_fitting_function(M, [5.46, -9.78, 4.17, -0.337, 0.57])
    end
end

"""
    Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95
"""
function KR13_acc(M::Float64)

    if M < 2.0
        return 0.0
    elseif 2.0 <= M <= 5.0
        param = [-0.0005950569221922047, 1.880258286365841e-5, 5.334076006529829 ]
        return (param[1] + param[2]*M^param[3])
    elseif 5.0 < M <= 15.0
        param = [-2.8696966498579606, 9.667563166507879, -8.877138312318019, 1.938386688261113, 0.1806112438315771]
        return kr_fitting_function(M, param)
    else
        return 0.21152
    end
end

"""
    Ryu et al. 2019, https://arxiv.org/abs/1905.04476
    values for 2.25 < M <= 5.0 extrapolated to entire range
"""
function Ryu19_acc(M::Float64)

    if M < 2.25
        return 0.0
    elseif M <= 34.0
        param = [-1.5255114554627316, 2.4026049650156693, -1.2534251472776456, 0.22152323784680614, 0.0335800899612107]
        return kr_fitting_function(M, param)
    else
        return 0.0348
    end

    # if 2.25 < M <= 5.0
    #     param = [1.1006004346467124, -2.923679036380424, 2.599624937871855, -0.9538130179325343, 0.16080793189704362]
    #     return (param[1] + param[2] * (M - 1.0) + param[3] * (M - 1.0)^2 + param[4] * (M - 1.0)^3 + param[5] * (M - 1.0)^4)/M^4
    # elseif M > 5.0
    #     param = [-2.8696966498579606, 9.667563166507879, -8.877138312318019, 1.938386688261113, 0.1806112438315771]
    #     return (param[1] + param[2] * (M - 1.0) + param[3] * (M - 1.0)^2 + param[4] * (M - 1.0)^3 + param[5] * (M - 1.0)^4)/M^4
    # else
    #     return 0.0
    # end
end

"""
    Caprioli&Spitkovsky 2015,
"""
function CS14_acc(M::Float64)
    vazza_factor = 0.5
    return vazza_factor * KR13_acc(M)
end

"""
    Constant efficiency as in Pfrommer+ 2016
"""
function P16_acc(M::Float64)
    return 0.5
end


"""
    fallback:
"""
function null_acc(M::Float64)
    return 0.0
end

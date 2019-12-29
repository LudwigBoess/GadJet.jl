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
        param = [-2.8696966498579606, 9.667563166507879,
                 -8.877138312318019, 1.938386688261113,
                 0.1806112438315771]
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
        param = [-1.5255114554627316, 2.4026049650156693,
                 -1.2534251472776456, 0.22152323784680614,
                  0.0335800899612107]
        return kr_fitting_function(M, param)
    else
        return 0.0348
    end
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

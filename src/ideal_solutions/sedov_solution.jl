using SpecialFunctions: gamma
using Interpolations: LinearInterpolation
using ProgressMeter

struct SedovData

     t::Float64
     m::Float64
     r::Vector{Float32}
     vr::Vector{Float32}
     rho::Vector{Float32}
     U::Vector{Float32}
     P::Vector{Float32}
     Pth::Vector{Float32}
     Pcr::Vector{Float32}
     E::Float64
     r_shock_rho::Float64
     r_shock_P::Float64
     rho_out::Float64
     rho_s::Float64
     cs_out::Float64
     gamma::Float64
     ndim::Int64

     function SedovData(t::Float64, m::Float64,
                        r::Vector{Float32}, vr::Vector{Float32},
                        rho::Vector{Float32}, U::Vector{Float32},
                        Pcr::Vector{Float32}, E::Float64;
                        γ::Float64=5.0/3.0, ndim::Int64=3)

        P = @. (γ - 1.0 ) * rho * U
        Pth = P

        if Pcr != zeros(Float32, length(P))
            P += Pcr
        end


        firstbin = Int(floor(length(rho) * 0.99)-1)

        rho_out = mean( rho[ firstbin:end ])
        P_out   = mean( P[ firstbin:end ] )

        cs_out  = sqrt( 5.0/3.0 * P_out/rho_out)

        rho_s = rho_out * ( γ + 1.0 )/( γ - 1.0 )

        k = findmax(rho)[2]
        r_shock_rho = r[k]

        k = findmax(P)[2]
        r_shock_P   = r[k]

        new(t, m, r, vr, rho, U, P, Pth, Pcr, E, r_shock_rho, r_shock_P,
            rho_out, rho_s, cs_out, γ, ndim)
    end

end

struct SedovFit

    alpha::Float64
    w::Int64
    D
    P
    V

end

struct SedovSolution

    r::Vector{Float32}
    vr::Vector{Float64}
    rho::Vector{Float64}
    P::Vector{Float64}
    U::Vector{Float64}
    r_shock::Float64
    v_shock::Float64
    xs::Float64
    mach::Float64

end

function R_s_analytic(i::Int64, data::SedovData, fit::SedovFit)
    return (data.E * data.t^2 / ( fit.alpha * data.rho[i] ))^(1.0/ ( 2.0 + data.ndim ))
end

function R_s_dot_analytic(i::Int64, data::SedovData, fit::SedovFit)
    t = maximum([data.t, 1.e-30])
    return 2.0 * R_s_analytic(i, data, fit) / ( (data.ndim + 2.0 - fit.w) * t)
end

function v_s_analytic(i::Int64, data::SedovData, fit::SedovFit)
    return 2.0 * R_s_dot_analytic(i, data, fit) /
            ( data.gamma + 1.0 )
end

function P_s_analytic(i::Int64, data::SedovData, fit::SedovFit)
    return 2.0 * data.rho_out * R_s_dot_analytic(i, data, fit)^2 / ( data.gamma + 1.0 )
end

function rho_analytic(i::Int64, data::SedovData, fit::SedovFit)
    R = R_s_analytic(i, data, fit)
    xi = data.r[i]/R
    if xi < 1.0
        return data.rho_s * fit.D(xi)
    else
        return data.rho_out
    end
end

function P_analytic(i::Int64, data::SedovData, fit::SedovFit)
    R = R_s_analytic(i, data, fit)
    xi = data.r[i]/R
    if xi < 1.0
        return P_s_analytic(i, data, fit) * fit.P(xi)
    else
        return 0.0
    end
end

function v_analytic(i::Int64, data::SedovData, fit::SedovFit)
    R = R_s_analytic(i, data, fit)
    xi = data.r[i]/R
    if xi < 1.0
        return v_s_analytic(i, data, fit) * fit.V(xi)
    else
        return 0.0
    end
end


struct fit_parameters

    ndim::Int64
    γ::Float64

    w::Float64
    w1::Float64
    w2::Float64
    w3::Float64

    b0::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    b4::Float64
    b5::Float64
    b6::Float64
    b7::Float64
    b8::Float64

    C0::Float64
    C1::Float64
    C2::Float64
    C3::Float64
    C4::Float64
    C5::Float64
    C6::Float64

end

struct w1Data
    ndim::Int64
    γ::Float64

    w::Float64
    w1::Float64
    w2::Float64
    w3::Float64

    b0::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    b4::Float64
    b5::Float64
    b6::Float64
    b7::Float64
    b8::Float64

    C0::Float64
    C1::Float64
    C2::Float64
    C3::Float64
    C4::Float64
    C5::Float64
    C6::Float64
end

struct w2Data
    ndim::Int64
    γ::Float64

    w::Float64
    w1::Float64
    w2::Float64
    w3::Float64

    b0::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    b4::Float64
    b5::Float64
    b6::Float64
    b7::Float64
    b8::Float64

    C0::Float64
    C1::Float64
    C2::Float64
    C3::Float64
    C4::Float64
    C5::Float64
    C6::Float64
end

struct w3Data
    ndim::Int64
    γ::Float64

    w::Float64
    w1::Float64
    w2::Float64
    w3::Float64

    b0::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    b4::Float64
    b5::Float64
    b6::Float64
    b7::Float64
    b8::Float64

    C0::Float64
    C1::Float64
    C2::Float64
    C3::Float64
    C4::Float64
    C5::Float64
    C6::Float64
end

struct elseData
    ndim::Int64
    γ::Float64

    w::Float64
    w1::Float64
    w2::Float64
    w3::Float64

    b0::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    b4::Float64
    b5::Float64
    b6::Float64
    b7::Float64
    b8::Float64

    C0::Float64
    C1::Float64
    C2::Float64
    C3::Float64
    C4::Float64
    C5::Float64
    C6::Float64
end

"""
        w1 functions
"""
function xifunc(F::Float64, c::w1Data)
    return F^(-c.b6) * (c.C1 * (c.F - c.C2))^c.b2 * (c.C3 * (c.C4 - F))^(-c.b1)
end
function Dfunc(F::Float64, c::w1Data)
    return xifunc(F, c)^(c.ndim-2)
end
function Vfunc(F::Float64, c::w1Data)
    return xifunc(F, c)
end
function Pfunc(F::Float64, c::w1Data)
    return xifunc(F, c)^c.ndim
end

"""
        w2 functions
"""
function xifunc(F::Float64, c::w2Data)
    return  F^(-c.b6) * (c.C1 * (F - c.C2))^((c.γ - 1.0 )*c.b0) *
            exp(((c.γ + 1.0)*c.b0*(1.0 - F)) / (F - c.C2))
end
function Dfunc(F::Float64, c::w2Data)
    return  F^c.b7 * (c.C1 * (F - c.C2))^((4.0 - c.ndim - 2.0*c.γ)*c.b0) *
            (c.C5 * (c.C6 - F))^(-c.b5) *
            exp((-2.0*(c.γ + 1.0)*c.b0*(1.0 - F)) / (F - c.C2))
end
function Vfunc(F::Float64, c::w2Data)
    return xifunc(F, c) * F
end
function Pfunc(F::Float64, c::w2Data)
    return F^c.b8 * (c.C1 * (F - c.C2))^(-c.ndim * c.γ * c.b0) *
           (c.C5 * (c.C6 - F))^(1.0 - c.b5)
end

"""
        w3 functions
"""
function xifunc(F::Float64, c::w3Data)
    return F^(-c.b6) * (c.C1 * (F - c.C2))^c.b2 * (c.C5 * (c.C6 - F))^(-c.b1)
end
function Dfunc(F::Float64, c::w3Data)
    return F^c.b7 * (c.C1 * (F - c.C2))^(c.b3-c.w*c.b2) *
           (c.C5 * (c.C6 - F))^(1.0 - 4.0*c.b0) *
           exp((-c.ndim*c.γ*(c.γ + 1.0)*c.b0*(1.0 - F)) /
           (c.C6 - F))
end
function Vfunc(F::Float64, c::w3Data)
    return xifunc(F, c) * F
end
function Pfunc(F::Float64, c::w3Data)
    return F^c.b8 * (c.C5 * (c.C6 - F))^(2.0*(c.ndim*c.γ-c.ndim-c.γ)*c.b0) *
           exp((-c.ndim*c.γ*(c.γ+1.0)*c.b0*(1.0 - F)) / (c.C6 - F))
end

"""
    else functions
"""
function xifunc(F::Float64, c::elseData)
    return F^(-c.b6) * (c.C1 * (F - c.C2))^c.b2 * (c.C3 * (c.C4 - F))^(-c.b1)
end
function Dfunc(F::Float64, c::elseData)
    return F^c.b7 * (c.C1 * (F - c.C2))^(c.b3 - c.w*c.b2) *
           (c.C3 * (c.C4 - F))^(c.b4 + c.w*c.b1) *
           (c.C5 * (c.C6 - F))^(-c.b5)
end
function Vfunc(F::Float64, c::elseData)
    return xifunc(F, c) * F
end
function Pfunc(F::Float64, c::elseData)
    return F^(c.b8) * (c.C3 * (c.C4 - F))^(c.b4 + (c.w - 2.0)*c.b1) *
           (c.C5 * (c.C6 - F))^(1.0 - c.b5)
end



function get_ideal_sedov_fit(rho_s::Float64; ndim::Int64=3, γ::Float64=5.0/3.0)

    w = 0 # w parametrizes the background density as rho = A * r^(-w)
    # The following is copied from Book (1991)
    # The omega parameters
    w1 = (3*ndim - 2 + γ*(2-ndim)) / (γ + 1)
    w2 = (2*(γ-1) + ndim) / γ
    w3 = ndim * (2-γ)
    # The beta_n parameters:
    b0 = 1 / (ndim * γ - ndim + 2)
    b2 = (γ - 1) / (γ*(w2-w))
    b3 = (ndim - w) / (γ*(w2-w))
    b5 = (2*ndim - w*(γ+1)) / (w3 - w)
    b6 = 2 / (ndim + 2 -w)
    b1 = b2 + (γ+1)*b0 -b6
    b4 = b1 * ((ndim - w) * (ndim + 2 - w)) / (w3 - w)
    b7 = w * b6
    b8 = ndim * b6
    # The C parameters
    C0 = 2^ndim * π^((ndim-1)/2) * gamma((ndim+1)/2) / gamma(ndim)
    C5 = 2 / (γ-1)
    C6 = (γ+1) / 2
    C1 = γ * C5
    C2 = C6 / γ
    C3 = (ndim*γ - ndim + 2) / ((w1-w)*C6)
    C4 = (ndim + 2 - w) * b0 * C6
    # Eq. (8) to (11) in Book (1991) or their later alternatives

    # println("w = $w\nw1 = $w1\nw2 = $w2")
    if w == w1
        # println("w1!")
        c = w1Data(ndim, γ,
                   w, w1, w2, w3,
                   b0, b1, b2, b3, b4, b5, b6, b7, b8,
                   C0, C1, C2, C3, C4, C5, C6)
    elseif w == w2
        # println("w2!")
        c = w2Data(ndim, γ,
                   w, w1, w2, w3,
                   b0, b1, b2, b3, b4, b5, b6, b7, b8,
                   C0, C1, C2, C3, C4, C5, C6)
    elseif w == w3
        # println("w3!")
        c = w3Data(ndim, γ,
                   w, w1, w2, w3,
                   b0, b1, b2, b3, b4, b5, b6, b7, b8,
                   C0, C1, C2, C3, C4, C5, C6)
    else
        # println("none!")
        c = elseData(ndim, γ,
                   w, w1, w2, w3,
                   b0, b1, b2, b3, b4, b5, b6, b7, b8,
                   C0, C1, C2, C3, C4, C5, C6)
    end


    if w1 > w
        Fmin = C2
    else
        Fmin = C6
    end

    # We use F as a proxy to do the integration over xi
    F = range(Fmin, stop=1.0, length=Int(1e5))
    xi = zeros(length(F))
    for i = 1:length(F)
        xi[i] = xifunc(F[i], c)
    end

    k = sortperm(F)
    F = F[k] # sort F according to the xi-values it procudes
    xi = xi[k]

    D = zeros(length(F))
    V = zeros(length(F))
    P = zeros(length(F))
    for i = 1:length(F)
        D[i] = Dfunc(F[i], c)
        V[i] = Vfunc(F[i], c)
        P[i] = Pfunc(F[i], c)
    end
    D_interpolated = LinearInterpolation(xi, D)
    P_interpolated = LinearInterpolation(xi, P)
    V_interpolated = LinearInterpolation(xi, V)
    I = @. xi^(ndim-1) * (D*V^2 + P) # Integrand
    dxi = zeros(length(I))
    # Use middle points
    for i = 1:(length(I)-1)
        I[i] = 0.5 * (I[i] + I[i+1])
        dxi[i] = xi[i+1] - xi[i]
    end
    integral = sum( @. I * dxi )
    # Finally, get alpha
    alpha = ((8 * C0) / ((γ^2 - 1)*(ndim + 2 - w)^2)) * integral

    return SedovFit(alpha, w, D_interpolated, P_interpolated, V_interpolated)
end

function get_ideal_sedov_data(data::SedovData, fit::SedovFit)

    N   = length(data.r)
    rho = zeros(N)
    P   = zeros(N)
    v   = zeros(N)

    for i = 1:N
        rho[i] = rho_analytic(i, data, fit)
        P[i]   = P_analytic(i, data, fit)
        v[i]   = v_analytic(i, data, fit)
    end

    U = @. P / ( (data.gamma - 1.0 ) * rho )

    k = findmax(rho)[2]

    r_shock = data.r[k]
    xs = rho[k]/rho[end]
    v_shock = v[k]#/(1.0 - 1.0/xs)
    mach = v_shock/data.cs_out

    return SedovSolution(data.r, v, rho, P, U, r_shock, v_shock, xs, mach)
end


function get_sedov_data_from_snapshot(fi::String, blast_center::Vector{Float64}=[3.0, 3.0, 3.0];
                                      CRs=false, Nbins::Int64=0)

    h = head_to_obj(fi)

    #info = read_info(fi)
    t = h.time
    m = h.massarr[1]

    x = read_block_by_name(fi, "POS", parttype=0)

    r = @. sqrt( (x[:,1] - blast_center[1])^2 +
                 (x[:,2] - blast_center[2])^2 +
                 (x[:,3] - blast_center[3])^2 )

    k = sortperm(r)

    r = Float32.(r[k])

    v = read_block_by_name(fi, "VEL", parttype=0)[k,:]

    vr = @. sqrt( v[:,1]^2 + v[:,2]^2 + v[:,3]^2 )

    rho = read_block_by_name(fi, "RHO", parttype=0)[k, 1]

    U = read_block_by_name(fi, "U", parttype=0)[k, 1]

    if CRs
        CRpP = read_block_by_name(fi, "CRpP", parttype=0)[k, 1]
        Ecr = @. CRpP/(1.0/3.0 * rho)
        γ = 7.0/5.0
    else
        CRpP = zeros(Float32, length(U))
        Ecr = zeros(Float32, length(U))
        γ = 5.0/3.0
    end

    E = sum( @. (Ecr * m + U * m + 0.5 * m * vr^2) )

    if Nbins != 0
        N = length(r)

        rbinwidth = maximum(r)/Nbins
        r_raw    = r
        v_raw    = vr
        rho_raw  = rho
        U_raw    = U
        CRpP_raw = CRpP

        r     = zeros(Float32, Nbins+1)
        vr    = zeros(Float32, Nbins+1)
        rho   = zeros(Float32, Nbins+1)
        U     = zeros(Float32, Nbins+1)
        CRpP  = zeros(Float32, Nbins+1)
        Npart = zeros(Float32, Nbins+1)

        @showprogress "Binning..." for i = 1:N
            bin = floor(Int64, r_raw[i]/rbinwidth) + 1

            vr[bin]     += v_raw[i]
            rho[bin]    += rho_raw[i]
            U[bin]      += U_raw[i]
            CRpP[bin]   += CRpP_raw[i]
            Npart[bin]  += 1
        end

        for i = 1:(Nbins+1)
            r[i] = (i-1) * rbinwidth
            vr[i]   /= Npart[i]
            rho[i]  /= Npart[i]
            U[i]    /= Npart[i]
            CRpP[i] /= Npart[i]
        end

    end

    k = findall(isnan.(rho))

    deleteat!(r, k)
    deleteat!(vr, k)
    deleteat!(rho, k)
    deleteat!(U, k)
    deleteat!(CRpP, k)


    return SedovData(t, m, r, vr, rho, U, CRpP, E, γ=γ)

end

function get_sedov_solution(filename::String, blast_center::Vector{Float64}=[3.0, 3.0, 3.0];
                            CRs::Bool=false, Nbins::Int64=0, Ndim::Int64=3)

    sedov_data = get_sedov_data_from_snapshot(filename, blast_center, CRs=CRs, Nbins=Nbins)

    if CRs
        γ=7.0/5.0
    else
        γ=5.0/3.0
    end

    sedov_fit  = get_ideal_sedov_fit(sedov_data.rho_s, ndim=Ndim, γ=γ)

    sedov_ideal = get_ideal_sedov_data(sedov_data, sedov_fit)

    return sedov_data, sedov_ideal

end

# Riemann solver with Cosmic Rays following Pfrommer et al 2016
# Use left and right state of initial condition to solve shocktube problems

using Roots


"""
    Helper functions to set up Machnumber dependent IC
"""
function solvePr(Pr::Float64,
                 rhol::Float64, rhor::Float64,
                 Pl::Float64, Mach::Float64,
                 γ::Float64)

    γ_1 = γ - 1.0
    γ_pow = γ_1/(2.0*γ)
    η2 = (γ-1.0)/(γ+1.0)

    c_l = sqrt.(γ*Pl/rhol)
    c_r = sqrt.(γ*Pr/rhor)

    f(x) = ( x/Pr - 1.0 ) * sqrt.( (1.0 - η2) / (γ * ( x/Pr + η2 ) )) -
            2.0 / γ_1 * c_l / c_r * ( 1.0 - ( x/Pl )^γ_pow )

    P_m = find_zero(f, (Pr, Pl), Bisection())


    ρ_mr = rhor * ( ( P_m + η2 * Pr ) / ( Pr + η2 * P_m))
    ρ_ml = rhol * ( P_m / Pl )^(1.0/γ)

    vm = 2.0 * c_l / γ_1 * ( 1.0 - (P_m/Pl)^γ_pow )

    vs = vm / ( 1.0 - rhor/ρ_mr)

    M = vs / c_r

    return M - Mach
end

function solvePrfromMach(rhol::Float64, rhor::Float64,
                         Pl::Float64, Mach::Float64,
                         γ::Float64)

    findPrHelper(P) = solvePr(P, rhol, rhor,
                              Pl, Mach, γ)

    Pr = find_zero(findPrHelper, (1.e-5, 1.e5), Bisection() )

    return Pr
end

function solveMach(Pl::Float64, Pr::Float64,
                   rhol::Float64, rhor::Float64,
                   γ::Float64)

   γ_1 = γ - 1.0
   γ_pow = γ_1/(2.0*γ)
   η2 = (γ-1.0)/(γ+1.0)

   c_l = sqrt.(γ*Pl/rhol)
   c_r = sqrt.(γ*Pr/rhor)

   f(x) = ( x/Pr - 1.0 ) * sqrt.( (1.0 - η2) / (γ * ( x/Pr + η2 ) )) -
           2.0 / γ_1 * c_l / c_r * ( 1.0 - ( x/Pl )^γ_pow )

   P_m = find_zero(f, (Pr, Pl), Bisection())


   ρ_mr = rhor * ( ( P_m + η2 * Pr ) / ( Pr + η2 * P_m))
   ρ_ml = rhol * ( P_m / Pl )^(1.0/γ)

   vm = 2.0 * c_l / γ_1 * ( 1.0 - (P_m/Pl)^γ_pow )

   vs = vm / ( 1.0 - rhor/ρ_mr)

   M = vs / c_r

   return M
end

"""
        DSA Models
"""

"""
    Kang&Ryu 2007, http://arxiv.org/abs/0704.1521v1
"""
function KR07(M::Float64)
    if M <= 2.0
        return 1.96e-3*(M^2 - 1.)               # eq. A3
    else
        b = [5.46, -9.78, 4.17, -0.337, 0.57]   # eq. A4
        η = 0.
        for i ∈ 1:length(b)
            η += b[i] * ((M - 1.)^(i-1))/M^4    # eq. A5
        end
        return η
    end
end

"""
    Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95
"""
function KR13(M::Float64)

    if M < 2.0
        return 0.0
    elseif 2.0 <= M <= 5.0
        param = [-0.0005950569221922047, 1.880258286365841e-5, 5.334076006529829 ]
        return param[1] + param[2]*M^param[3]
    elseif 5.0 < M <= 15.0
        param = [-2.8696966498579606, 9.667563166507879, -8.877138312318019, 1.938386688261113, 0.1806112438315771]
        return (param[1] + param[2] * (M - 1.0) + param[3] * (M - 1.0)^2 + param[4] * (M - 1.0)^3 + param[5] * (M - 1.0)^4)/M^4
    else
        return 0.21152
    end
end

"""
    Ryu et al. 2019, https://arxiv.org/abs/1905.04476
    values for 2.25 < M <= 5.0 extrapolated to entire range
"""
function R_19(M::Float64)

    if M < 2.25
        return 0.0
    elseif M <= 34.0
        param = [-1.5255114554627316, 2.4026049650156693, -1.2534251472776456, 0.22152323784680614, 0.0335800899612107]
        return (param[1] + param[2] * (M - 1.0) + param[3] * (M - 1.0)^2 + param[4] * (M - 1.0)^3 + param[5] * (M - 1.0)^4)/M^4
    else
        return  0.0348
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
function CS15(M::Float64)
    vazza_factor = 0.5
    return vazza_factor * KR13(M)
end

function null_eff(M::Float64)
    return 0.0
end

"""
    Datatypes for IC Parameters and Solution
"""
mutable struct RiemannParameters#{F<:Function}
    rhol::Float64
    rhor::Float64
    Pl::Float64
    Pr::Float64
    Ul::Float64
    Ur::Float64
    cl::Float64
    cr::Float64
    M::Float64
    t::Float64
    x_contact::Float64
    Pe_ratio::Float64
    γ_th::Float64
    γ_cr::Float64
    γ_exp::Float64
    η2::Float64
    eff_function#::F

    function RiemannParameters(;rhol::Float64=1.0,  rhor::Float64=0.125,
                                Pl::Float64=0.0,    Pr::Float64=0.0,
                                Ul::Float64=0.0,    Ur::Float64=0.0,
                                Mach::Float64=0.0,  t::Float64,
                                x_contact::Float64=70.0,
                                Pe_ratio::Float64=0.01,
                                γ_th::Float64=5.0/3.0,
                                γ_cr::Float64=4.0/3.0,
                                eff_model::Int64=0)

        γ_exp    = ( γ_th - 1.0 )/( 2.0 * γ_th )
        η2       = (γ_th-1.0)/(γ_th+1.0)

        # calculate Ul and Pl depending on input
        if (Pl == 0.0) & (Ul != 0.0)
            Pl = ( γ_th - 1.0 ) * rhol * Ul
        elseif (Ul == 0.0) & (Pl != 0.0)
            Ul = Pl / ( (γ_th - 1.0) * rhol )
        else
            println("Error! Both Ul and Pl are zero!")
        end

        # calculate Ur and Pr depending on input
        if (Pr == 0.0) & (Ur != 0.0)
            Pr = ( γ_th - 1.0 ) * rhor * Ur
        elseif (Ur == 0.0) & (Pr != 0.0)
            Ur = Pr / ( (γ_th - 1.0) * rhor )
        elseif (Ur == 0.0) & (Pr == 0.0) & (Mach == 0.0)
            println("Error! Ur, Pr and Mach are zero! Can't find solution!")
        else
            println("Both Ur and Pr are zero! Will calculate them depending on Machnumber.")
            Pr = solvePrfromMach(rhol, rhor, Pl, Mach, γ_th)
            Ur = Pr / ( (γ_th - 1.0) * rhor )
        end

        if Mach == 0.0
            Mach = solveMach(Pl, Pr, rhol, rhor, γ_th)
        end

        cl = sqrt.(γ_th * Pl / rhol)
        cr = sqrt.(γ_th * Pr / rhor)

        if eff_model == 0
            eff_function = null_eff
        elseif eff_model == 1
            eff_function = KR07
        elseif eff_model == 2
            eff_function = KR13
        elseif eff_model == 3
            eff_function = R_19
        elseif eff_model == 4
            eff_function = CS15
        else
            println("Invalid DSA model selection!")
            println("Selecting null function")
            eff_function = null_eff
        end

        #new{typeof(F)}
        new(rhol, rhor,
            Pl, Pr,
            Ul, Ur,
            cl, cr,
            Mach, t,
            x_contact,
            Pe_ratio,
            γ_th, γ_cr, γ_exp, η2,
            eff_function)
    end
end

mutable struct RiemannSolution
    x::Array{Float64,1}
    rho::Array{Float64,1}
    rho4::Float64
    rho3::Float64
    P_tot::Array{Float64,1}
    P_th::Array{Float64,1}
    P_cr::Array{Float64,1}
    P_cr_p::Array{Float64,1}
    P_cr_e::Array{Float64,1}
    P34::Float64
    U_tot::Array{Float64,1}
    U_th::Array{Float64,1}
    E_cr_p::Array{Float64,1}
    E_cr_e::Array{Float64,1}
    v::Array{Float64,1}
    v34::Float64
    vt::Float64
    vs::Float64
    Mach::Float64

    function RiemannSolution(x::Array{Float64,1})
        N = length(x)
        new(x,
            zeros(N),   # rho
            0.0,        # rho2
            0.0,        # rho3
            zeros(N),   # P_tot
            zeros(N),   # P_th
            zeros(N),   # P_cr
            zeros(N),   # P_cr_e
            zeros(N),   # P_cr_p
            0.0,        # P34
            zeros(N),   # U_tot
            zeros(N),   # U_th
            zeros(N),   # E_cr_p
            zeros(N),   # E_cr_e
            zeros(N),   # v
            0.0,        # v23
            0.0,        # vt
            0.0,        # vs
            0.0)        # Mach
    end
end



"""
    Functions to solve individual segments of the shocktube
"""


"""
    Pressure
"""
function solveP34(;par::RiemannParameters)
    # Solves the pressure of the middle state in a shocktube

    P_m(P) = ( P/par.Pr - 1.0 ) * sqrt.( ( 1.0 - par.η2 ) / ( par.γ_th *
             ( P/par.Pr + par.η2 ) )) - 2.0 / (par.γ_th - 1.0 ) * par.cl / par.cr *
             ( 1.0 - ( P/par.Pl )^par.γ_exp )

    P34 = find_zero( P_m, (par.Pr, par.Pl), Bisection() )

    return P34
end

function solveP2(x::Float64; par::RiemannParameters)
    # solves the pressure along the rarefaction wave
    return par.Pl * ( -par.η2 * x / ( par.cl * par.t ) + ( 1.0 - par.η2 ) )^( 1.0/par.γ_exp )
end

function solveP(x::Float64; par::RiemannParameters, sol::RiemannSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return par.Pl, par.Pl, 0.0, 0.0, 0.0
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        P2 = solveP2(x, par=par)
        return P2, P2, 0.0, 0.0, 0.0
    elseif -sol.vt * par.t < x <= sol.v34 * par.t
        return sol.P34, sol.P34, 0.0, 0.0, 0.0
    elseif sol.v34 * par.t < x <= sol.vs * par.t
        η_cr = par.eff_function(par.M)
        return sol.P34,                             # total pressure
               (1.0 - η_cr)*sol.P34,                # thermal pressure
               η_cr*sol.P34,                        # total cr pressure
               (1.0 - par.Pe_ratio)*η_cr*sol.P34,   # cr proton pressure
               par.Pe_ratio*η_cr*sol.P34            # cr electron pressure
    elseif sol.vs * par.t < x
        return par.Pr, par.Pr, 0.0, 0.0, 0.0
    else
        println("Error!")
        return 0.0, 0.0, 0.0, 0.0, 0.0
    end
end


"""
    Density
"""
function solveRho2(x::Float64 ;par::RiemannParameters)
    # solves density along the rarefaction wave
    return par.rhol * ( -par.η2 * x / ( par.cl * par.t ) +
                      ( 1.0 - par.η2 ) )^( 2.0 / (par.γ_th - 1.0 ) )
end

function solveRho3(;par::RiemannParameters, sol::RiemannSolution)
    # solves density left of contact discontiuity and right of refraction wave
    return par.rhol * ( sol.P34/par.Pl )^( 1.0/par.γ_th )
end

# function solveRho4(;par::RiemannParameters)
#     # solves density left of shockfront
#     return ( ( par.γ_th + 1.0 ) * par.M^2 ) /
#            ( 2.0 + ( par.γ_th - 1.0 ) * par.M^2 ) * par.rhor
# end

function solveRho4(;par::RiemannParameters, sol::RiemannSolution)

    η_cr = par.eff_function(par.M)
    P_cr = η_cr*sol.P34

    rho0 = ( ( par.γ_th + 1.0 ) * par.M^2 ) /
           ( 2.0 + ( par.γ_th - 1.0 ) * par.M^2 ) * par.rhor

    C = 2.0*(par.γ_th - par.γ_cr)*P_cr/((par.γ_th - 1.0)*sol.v34^2*(par.γ_cr - 1.0))

    find_rho1(rho1) = rho0/(1.0 - C/rho1) - rho1

    ρ1 = find_zero( find_rho1, (0.9*rho0, 10.0*rho0), Bisection() )

    return ρ1

end

function solveRho(x::Float64; par::RiemannParameters, sol::RiemannSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return par.rhol
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return solveRho2(x, par=par)
    elseif -sol.vt * par.t < x <= sol.v34 * par.t
        return sol.rho3
    elseif sol.v34 * par.t < x <= sol.vs * par.t
        return sol.rho4
    elseif sol.vs * par.t < x
        return par.rhor
    else
        println("Error!")
        return 0.0
    end
end


"""
    Velocity
"""
function solveV34(;par::RiemannParameters, sol::RiemannSolution)
    # solves velocity along the isopressure region 2-3
    return 2.0 * par.cl / ( par.γ_th - 1.0 ) * ( 1.0 - (sol.P34 / par.Pl)^par.γ_exp )
end

function solveV2(x::Float64; par::RiemannParameters)
    # solves the velocity along the rarefaction wave
    return ( 1.0 - par.η2 ) * ( x / par.t + par.cl )
end

function solveVs(;par::RiemannParameters, sol::RiemannSolution)
    # solves the shock velocity
    return sol.v34 / ( 1.0 - par.rhor / sol.rho4 )
end

function solveVt(;par::RiemannParameters, sol::RiemannSolution)
    # solves the velocity of the tail of the rarefaction wave
    return par.cl - sol.v34/( 1.0 - par.η2 )
end

function solveV(x::Float64; par::RiemannParameters, sol::RiemannSolution)
    # solves the velocity along the shocktube
    if x <= -par.cl * par.t
        return 0.0
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return solveV2(x, par=par)
    elseif -sol.vt * par.t < x <= sol.vs * par.t
        return sol.v34
    elseif sol.vs * par.t < x
        return 0.0
    else
        println("Error!")
        return 0.0
    end
end


"""
    Main function for solver
"""
function solveHydroShock(x::Array{Float64,1}; par::RiemannParameters)
    # solves a standard sod shock with CR acceleration

    # set up datatype to store riemann solution
    sol = RiemannSolution(x)

    # transform into rest-frame of contact discontiuity
    x_in = sol.x .- par.x_contact

    # solve Pressure
    sol.P34 = solveP34(par=par)

    # solve velocity
    sol.v34 = solveV34(par=par, sol=sol)

    sol.rho4 = solveRho4(par=par, sol=sol)

    sol.vs = solveVs(par=par, sol=sol)
    sol.vt = solveVt(par=par, sol=sol)
    sol.v = solveV.(x_in, par=par, sol=sol)

    for i = 1:length(sol.x)
        sol.P_tot[i], sol.P_th[i], sol.P_cr[i], sol.P_cr_p[i], sol.P_cr_e[i] = solveP(x_in[i], par=par, sol=sol)
    end

    # solve density
    sol.rho3 = solveRho3(par=par, sol=sol)

    for i = 1:length(sol.x)
        sol.rho[i] = solveRho(x_in[i], par=par, sol=sol)
    end

    sol.Mach = par.M

    sol.U_tot = sol.P_tot ./ ( (par.γ_th - 1.0) .* sol.rho )
    sol.U_th = sol.P_th ./ ( (par.γ_th - 1.0) .* sol.rho )
    sol.E_cr_p = sol.P_cr_p ./ ( (par.γ_cr - 1.0) .* sol.rho )
    sol.E_cr_e = sol.P_cr_e ./ ( (par.γ_cr - 1.0) .* sol.rho )

    return sol
end

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
    Datatypes for IC Parameters and Solution
"""
mutable struct RiemannParameters
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

    function RiemannParameters(;rhol::Float64=1.0, rhor::Float64=0.125,
                                Pl::Float64=0.0, Pr::Float64=0.0,
                                Ul::Float64=0.0, Ur::Float64=0.0,
                                Mach::Float64=0.0, t::Float64,
                                x_contact::Float64=70.0,
                                Pe_ratio::Float64=0.01,
                                γ_th::Float64=5.0/3.0,
                                γ_cr::Float64=4.0/3.0)


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

        new(rhol, rhor,
            Pl, Pr,
            Ul, Ur,
            cl, cr,
            Mach, t,
            x_contact,
            Pe_ratio,
            γ_th, γ_cr, γ_exp, η2)
    end
end

mutable struct RiemannSolution
    x::Array{Float64,1}
    rho::Array{Float64,1}
    rho2::Float64
    rho3::Float64
    P_tot::Array{Float64,1}
    P_th::Array{Float64,1}
    P_cr::Array{Float64,1}
    P_cr_p::Array{Float64,1}
    P_cr_e::Array{Float64,1}
    P23::Float64
    U::Array{Float64,1}
    v::Array{Float64,1}
    v23::Float64
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
            0.0,        # P23
            zeros(N),   # U
            zeros(N),   # v
            0.0,        # v23
            0.0,        # vt
            0.0,        # vs
            0.0)        # Mach
    end
end

"""
    Acceleration efficiencies for CRs
        (Kang et al. 2006, http://arxiv.org/abs/0704.1521v1)
"""
function ηle2(M)
    return 1.96e-3*(M^2 - 1.)               # eq. A3
end

function ηg2(M)
    b = [5.46, -9.78, 4.17, -0.337, 0.57]   # eq. A4
    η = 0.
    for i ∈ 1:length(b)
        η += b[i] * ((M - 1.)^(i-1))/M^4    # eq. A5
    end
    return η
end

function get_eff(M)
    if M <= 2.0
        return ηle2(M)
    else
        return ηg2(M)
    end
end


"""
    Functions to solve individual segments of the shocktube
"""


"""
    Pressure
"""
function solveP23(;par::RiemannParameters)
    # Solves the pressure of the middle state in a shocktube

    P_m(P) = ( P/par.Pr - 1.0 ) * sqrt.( ( 1.0 - par.η2 ) / ( par.γ_th *
             ( P/par.Pr + par.η2 ) )) - 2.0 / (par.γ_th - 1.0 ) * par.cl / par.cr *
             ( 1.0 - ( P/par.Pl )^par.γ_exp )

    P23 = find_zero( P_m, (par.Pr, par.Pl), Bisection() )

    return P23
end

function solveP4(x::Float64; par::RiemannParameters)
    # solves the pressure along the rarefaction wave
    return par.Pl * ( -par.η2 * x / ( par.cl * par.t ) + ( 1.0 - par.η2 ) )^( 1.0/par.γ_exp )
end

function solveP(x::Float64; par::RiemannParameters, sol::RiemannSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return par.Pl, par.Pl, 0.0, 0.0, 0.0
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        P4 = solveP4(x, par=par)
        return P4, P4, 0.0, 0.0, 0.0
    elseif -sol.vt * par.t < x <= sol.v23 * par.t
        return sol.P23, sol.P23, 0.0, 0.0, 0.0
    elseif sol.v23 * par.t < x <= sol.vs * par.t
        η_cr = get_eff(par.M)
        return sol.P23,                             # total pressure
               (1.0 - η_cr)*sol.P23,                # thermal pressure
               η_cr*sol.P23,                        # total cr pressure
               (1.0 - par.Pe_ratio)*η_cr*sol.P23,   # cr proton pressure
               par.Pe_ratio*η_cr*sol.P23            # cr electron pressure
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
function solveRho4(x::Float64 ;par::RiemannParameters)
    # solves density along the rarefaction wave
    return par.rhol * ( -par.η2 * x / ( par.cl * par.t ) +
                      ( 1.0 - par.η2 ) )^( 2.0 / (par.γ_th - 1.0 ) )
end

function solveRho3(;par::RiemannParameters, sol::RiemannSolution)
    # solves density left of contact discontiuity and right of refraction wave
    return par.rhol * ( sol.P23/par.Pl )^( 1.0/par.γ_th )
end

function solveRho2(;par::RiemannParameters)
    # solves density left of shockfront
    return ( ( par.γ_th + 1.0 ) * par.M^2 ) /
           ( 2.0 + ( par.γ_th - 1.0 ) * par.M^2 ) * par.rhor
end

function solveRho(x::Float64; par::RiemannParameters, sol::RiemannSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return par.rhol
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return solveRho4(x, par=par)
    elseif -sol.vt * par.t < x <= sol.v23 * par.t
        return sol.rho3
    elseif sol.v23 * par.t < x <= sol.vs * par.t
        return sol.rho2
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
function solveV23(;par::RiemannParameters, sol::RiemannSolution)
    # solves velocity along the isopressure region 2-3
    return 2.0 * par.cl / ( par.γ_th - 1.0 ) * ( 1.0 - (sol.P23 / par.Pl)^par.γ_exp )
end

function solveV4(x::Float64; par::RiemannParameters)
    # solves the velocity along the rarefaction wave
    return ( 1.0 - par.η2 ) * ( x / par.t + par.cl )
end

function solveVs(;par::RiemannParameters, sol::RiemannSolution)
    # solves the shock velocity
    return sol.v23 / ( 1.0 - par.rhor / sol.rho2 )
end

function solveVt(;par::RiemannParameters, sol::RiemannSolution)
    # solves the velocity of the tail of the rarefaction wave
    return par.cl - sol.v23/( 1.0 - par.η2 )
end

function solveV(x::Float64; par::RiemannParameters, sol::RiemannSolution)
    # solves the velocity along the shocktube
    if x <= -par.cl * par.t
        return 0.0
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return solveV4(x, par=par)
    elseif -sol.vt * par.t < x <= sol.vs * par.t
        return sol.v23
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
    sol.P23 = solveP23(par=par)

    # solve velocity
    sol.v23 = solveV23(par=par, sol=sol)

    sol.rho2 = solveRho2(par=par)

    sol.vs = solveVs(par=par, sol=sol)
    sol.vt = solveVt(par=par, sol=sol)
    sol.v = solveV.(x_in, par=par, sol=sol)

    for i = 1:length(sol.x)
        sol.P_tot[i], sol.P_th[i], sol.P_cr[i], sol.P_cr_p[i], sol.P_cr_e[i] = solveP(x_in[i], par=par, sol=sol)
    end

    # solve density

    sol.rho3 = solveRho3(par=par, sol=sol)
    sol.rho = solveRho.(x_in, par=par, sol=sol)

    sol.Mach = par.M

    sol.U = sol.P_tot ./ ( (par.γ_th - 1.0) .* sol.rho )

    return sol
end

# Riemann solver for standard hydro shocktube following Pfrommer et al 2006
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
    #ρ_ml = rhol * ( P_m / Pl )^(1.0/γ)

    vm = 2.0 * c_l / γ_1 * ( 1.0 - (P_m/Pl)^γ_pow )

    vs = vm / ( 1.0 - rhor/ρ_mr)

    M = vs / c_r

    return M - Mach
end

function solvePrfromMach(rhol::Float64, rhor::Float64,
                         Pl::Float64, Mach::Float64,
                         γ::Float64)

    # first approx without CRs
    findPrHelper(P) = solvePr(P, rhol, rhor,
                              Pl, Mach, γ)

    Pr = find_zero(findPrHelper, (1.e-5, 1.e5), Bisection() )


    return Pr
end

function solvePlfromMach(rhol::Float64, rhor::Float64,
                         Pr::Float64, Mach::Float64,
                         γ::Float64)

    # first approx without CRs
    findPlHelper(P) = solvePr(Pr, rhol, rhor,
                              P, Mach, γ)

    Pl = find_zero(findPlHelper, (1.e-5, 1.e5), Bisection() )


    return Pl
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
mutable struct SodParameters

    rhol::Float64           # denisty left
    rhor::Float64           # density right
    Pl::Float64             # pressure left
    Pr::Float64             # pressure right
    Ul::Float64             # internal energy left
    Ur::Float64             # internal energy right
    cl::Float64             # soundspeed left
    cr::Float64             # soundspeed right
    M::Float64              # Mach number
    t::Float64              # time
    x_contact::Float64      # position of the contact discontinuity along the tube
    γ_th::Float64           # adiabatic index of the gas
    γ_exp::Float64          # helper variable
    η2::Float64             # helper variable

    function SodParameters(;rhol::Float64=1.0,  rhor::Float64=0.125,
                            Pl::Float64=0.0,    Pr::Float64=0.0,
                            Ul::Float64=0.0,    Ur::Float64=0.0,
                            Mach::Float64=0.0,  t::Float64,
                            x_contact::Float64=70.0,
                            γ_th::Float64=5.0/3.0)

        γ_exp    = ( γ_th - 1.0 )/( 2.0 * γ_th )
        η2       = (γ_th-1.0)/(γ_th+1.0)

        # calculate Ul and Pl depending on input
        if (Pl == 0.0) & (Ul != 0.0)
            Pl = ( γ_th - 1.0 ) * rhol * Ul
        elseif (Ul == 0.0) & (Pl != 0.0)
            Ul = Pl / ( (γ_th - 1.0) * rhol )
        end

        # calculate Ur and Pr depending on input
        if (Pr == 0.0) & (Ur != 0.0)
            Pr = ( γ_th - 1.0 ) * rhor * Ur
        elseif (Ur == 0.0) & (Pr != 0.0)
            Ur = Pr / ( (γ_th - 1.0) * rhor )
        end

        # error handling
        if (Pr == 0.0) & (Pl == 0.0)
            error("No initial Pressure or energy values given!")
        end

        # solve right or left initial state from target mach number
        if (Pr == 0.0) & (Mach != 0.0)
            Pr = solvePrfromMach(rhol, rhor, Pl, Mach, γ_th)
            Ur = Pr / ( (γ_th - 1.0) * rhor )
        elseif (Pl == 0.0) & (Mach != 0.0)
            Pl = solvePlfromMach(rhol, rhor, Pr, Mach, γ_th)
            Ul = Pl / ( (γ_th - 1.0) * rhol )
        elseif (Ur == 0.0) & (Pr == 0.0) & (Mach == 0.0)
            error("Ur, Pr and Mach are zero! Can't find solution!")
        elseif (Ul == 0.0) & (Pl == 0.0) & (Mach == 0.0)
            error("Ul, Pl and Mach are zero! Can't find solution!")
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
            γ_th, γ_exp, η2)
    end
end

mutable struct SodHydroSolution
    x::Array{Float64,1}         # array of given positions
    rho::Array{Float64,1}       # array of densities along the tube
    rho4::Float64               # density in postshock region
    rho3::Float64               # density between contact disc. and rarefaction wave
    P::Array{Float64,1}         # array of pressures along the tube
    P34::Float64                # pressure between shock and rarefaction wave
    U::Array{Float64,1}         # array of internal energies along the tube
    v::Array{Float64,1}         # array of velocities along the tube
    v34::Float64                # velocity between shock and rarefaction wave
    vt::Float64                 # velocity of rarefaction wave
    vs::Float64                 # shock velocity
    Mach::Float64               # Mach number

    function SodHydroSolution(x::Array{Float64,1})
        N = length(x)
        new(x,
            zeros(N),   # rho
            0.0,        # rho2
            0.0,        # rho3
            zeros(N),   # P
            0.0,        # P34
            zeros(N),   # U
            zeros(N),   # v
            0.0,        # v23
            0.0,        # vt
            0.0,        # vs
            0.0)        # Mach
    end

    # multiple dispatch, in case gadget data is passed directly to the constructor
    # function SodHydroSolution(x::Array{Float32,1})
    #     N = length(x)
    #     new(x,
    #         zeros(N),   # rho
    #         0.0,        # rho2
    #         0.0,        # rho3
    #         zeros(N),   # P
    #         0.0,        # P34
    #         zeros(N),   # U
    #         zeros(N),   # v
    #         0.0,        # v23
    #         0.0,        # vt
    #         0.0,        # vs
    #         0.0)        # Mach
    # end
end



"""
    Functions to solve individual segments of the shocktube
"""


"""
    Pressure
"""
function solveP34(par::SodParameters)
    # Solves the pressure of the middle state in a shocktube

    P_m(P) = ( P/par.Pr - 1.0 ) * sqrt.( ( 1.0 - par.η2 ) / ( par.γ_th *
             ( P/par.Pr + par.η2 ) )) - 2.0 / (par.γ_th - 1.0 ) * par.cl / par.cr *
             ( 1.0 - ( P/par.Pl )^par.γ_exp )

    P34 = find_zero( P_m, (par.Pr, par.Pl), Bisection() )

    return P34
end

function solveP2(x::Float64, par::SodParameters)
    # solves the pressure along the rarefaction wave
    return par.Pl * ( -par.η2 * x / ( par.cl * par.t ) + ( 1.0 - par.η2 ) )^( 1.0/par.γ_exp )
end

function solveP(x::Float64, par::SodParameters, sol::SodHydroSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return par.Pl
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        P2 = solveP2(x, par)
        return P2
    elseif -sol.vt * par.t < x <= sol.v34 * par.t
        return sol.P34
    elseif sol.v34 * par.t < x <= sol.vs * par.t
        return sol.P34
    elseif sol.vs * par.t < x
        return par.Pr
    else
        println("Error!")
        return 0.0
    end
end


"""
    Density
"""
function solveRho2(x::Float64, par::SodParameters)
    # solves density along the rarefaction wave
    return par.rhol * ( -par.η2 * x / ( par.cl * par.t ) +
                      ( 1.0 - par.η2 ) )^( 2.0 / (par.γ_th - 1.0 ) )
end

function solveRho3(par::SodParameters, sol::SodHydroSolution)
    # solves density left of contact discontiuity and right of refraction wave
    return par.rhol * ( sol.P34/par.Pl )^( 1.0/par.γ_th )
end

function solveRho4(par::SodParameters)
    # solves density left of shockfront
    return ( ( par.γ_th + 1.0 ) * par.M^2 ) /
           ( 2.0 + ( par.γ_th - 1.0 ) * par.M^2 ) * par.rhor
end


function solveRho(x::Float64, par::SodParameters, sol::SodHydroSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return par.rhol
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return solveRho2(x, par)
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
function solveV34(par::SodParameters, sol::SodHydroSolution)
    # solves velocity along the isopressure region 2-3
    return 2.0 * par.cl / ( par.γ_th - 1.0 ) * ( 1.0 - (sol.P34 / par.Pl)^par.γ_exp )
end

function solveV2(x::Float64, par::SodParameters)
    # solves the velocity along the rarefaction wave
    return ( 1.0 - par.η2 ) * ( x / par.t + par.cl )
end

function solveVs(par::SodParameters, sol::SodHydroSolution)
    # solves the shock velocity
    return sol.v34 / ( 1.0 - par.rhor / sol.rho4 )
end

function solveVt(par::SodParameters, sol::SodHydroSolution)
    # solves the velocity of the tail of the rarefaction wave
    return par.cl - sol.v34/( 1.0 - par.η2 )
end

function solveV(x::Float64, par::SodParameters, sol::SodHydroSolution)
    # solves the velocity along the shocktube
    if x <= -par.cl * par.t
        return 0.0
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return solveV2(x, par)
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
function solveSodShock(x::Array{Float64,1}; par::SodParameters)
    # solves a standard sod shock with zero initial velocity

    N = length(x)

    # set up datatype to store riemann solution
    sol = SodHydroSolution(x)

    # transform into rest-frame of contact discontiuity
    x_in = sol.x .- par.x_contact

    # solve Pressure
    sol.P34 = solveP34(par)

    sol.rho3 = solveRho3(par, sol)
    sol.rho4 = solveRho4(par)

    # solve velocity
    sol.v34 = solveV34(par, sol)
    sol.vs = solveVs(par, sol)
    sol.vt = solveVt(par, sol)
    for i = 1:N
        sol.v[i] = solveV(x_in[i], par, sol)
    end

    for i = 1:N
        sol.P[i] = solveP(x_in[i], par, sol)
    end

    # solve density


    for i = 1:N
        sol.rho[i] = solveRho(x_in[i], par, sol)
    end

    sol.Mach = par.M

    sol.U = sol.P ./ ( (par.γ_th - 1.0) .* sol.rho )

    return sol
end

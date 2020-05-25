#include("cr_sod_shock_main.jl")

using Roots

mutable struct SodCRParameters_noCRs

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
    acc_function
    ξ::Float64
    first_guess::Float64

    function SodCRParameters_noCRs(;rhol::Float64=1.0,  rhor::Float64=0.125,
                                Pl::Float64=0.0,    Pr::Float64=0.0,
                                Ul::Float64=0.0,    Ur::Float64=0.0,
                                Mach::Float64=0.0,  t::Float64,
                                x_contact::Float64=70.0,
                                Pe_ratio::Float64=0.01,
                                γ_th::Float64=5.0/3.0,
                                γ_cr::Float64=4.0/3.0,
                                thetaB::Float64=0.0,
                                theta_crit::Float64=(π/4.0),
                                dsa_model::Int64=-1,
                                xs_first_guess::Float64=4.7)

        γ_exp    = ( γ_th - 1.0 )/( 2.0 * γ_th )
        η2       = (γ_th-1.0)/(γ_th+1.0)

        # calculate Ul and Pl depending on input
        if (Pl == 0.0) & (Ul != 0.0)
            Pl = ( γ_th - 1.0 ) * rhol * Ul
        elseif (Ul == 0.0) & (Pl != 0.0)
            Ul = Pl / ( (γ_th - 1.0) * rhol )
        else
            error("Both Ul and Pl are zero!")
        end


        if dsa_model == 0
            acc_function = KR07_acc
        elseif dsa_model == 1
            acc_function = KR13_acc
        elseif dsa_model == 2
            acc_function = Ryu19_acc
        elseif dsa_model == 3
            acc_function = CS14_acc
        elseif dsa_model == 4
            acc_function = P16_acc
        else
            error("Invalid DSA model selection!")

        end

        # calculate Ur and Pr depending on input
        if (Pr == 0.0) & (Ur != 0.0)
            Pr = ( γ_th - 1.0 ) * rhor * Ur
        elseif (Ur == 0.0) & (Pr != 0.0)
            Ur = Pr / ( (γ_th - 1.0) * rhor )
        elseif (Ur == 0.0) & (Pr == 0.0) & (Mach == 0.0)
            error("Ur, Pr and Mach are zero! Can't find solution!")
        else
            println("Both Ur and Pr are zero! Will calculate them depending on Machnumber.")
            Pr = solvePrfromMachCR(Pl, Pr,
                                   rhor, rhol,
                                   γ_th, γ_cr, γ_exp,
                                   Mach, acc_function)
            Ur = Pr / ( (γ_th - 1.0) * rhor )
        end

        # calculate B angle dependent efficiency following Pais+ 2018, MNRAS, 478, 5278
        delta_theta = π/18.0
        thetaB *= (π/180.0)
        etaB = 0.5*( tanh( (theta_crit - thetaB)/delta_theta ) + 1.0 )

        ξ = etaB*acc_function(Mach)/(1.0 - etaB*acc_function(Mach))

        # if dsa_model == 4
        #     ξ = etaB*acc_function(Mach)/(1.0 - etaB*acc_function(Mach))
        # else
        #     ξ = etaB*acc_function(Mach)
        # end
        #ξ = acc_function(Mach)/(1.0 - acc_function(Mach))


        cl = sqrt.( γ_th * Pl / rhol)
        cr = sqrt.( γ_th * Pr / rhor)


        new(rhol, rhor,
            Pl, Pr,
            Ul, Ur,
            cl, cr,
            Mach, t,
            x_contact,
            Pe_ratio,
            γ_th, γ_cr,
            γ_exp,
            η2,
            acc_function,
            ξ,
            xs_first_guess)

    end
end


# Riemann solver for standard hydro shocktube following Pfrommer et al 2006
# Use left and right state of initial condition to solve shocktube problems

"""
    Helper functions to set up Machnumber dependent IC
"""
function solvePrCR(Pr::Float64,
                 rhol::Float64, rhor::Float64,
                 Pl::Float64, Mach::Float64,
                 γ::Float64, γ_cr::Float64,
                 eff_function)

    γ_1 = γ - 1.0
    γ_pow = γ_1/(2.0*γ)
    η2 = (γ-1.0)/(γ+1.0)

    c_l = sqrt.(γ*Pl/rhol)
    c_r = sqrt.(γ*Pr/rhor)

    f(x) = ( x/Pr - 1.0 ) * sqrt.( (1.0 - η2) / (γ * ( x/Pr + η2 ) )) -
            2.0 / γ_1 * c_l / c_r * ( 1.0 - ( x/Pl )^γ_pow )

    P_m = find_zero(f, (Pr, Pl), Bisection())

    # CR part

    P_cr = eff_function(Mach) * P_m

    vm = 2.0 * c_l / γ_1 * ( 1.0 - (P_m/Pl)^γ_pow )

    C = 2.0*(γ - γ_cr)*P_cr/( γ_1* vm^2 * (γ_cr - 1.0) )

    ρ_mr0 = rhor * ( ( P_m + η2 * Pr ) / ( Pr + η2 * P_m))
    #ρ_ml = rhol * ( P_m / Pl )^(1.0/γ)

    find_rho1(rho1) = ρ_mr0/(1.0 - C/rho1) - rho1

    ρ_mr = find_zero( find_rho1, ρ_mr0 )

    vm = 2.0 * c_l / γ_1 * ( 1.0 - (P_m/Pl)^γ_pow )

    vs = vm / ( 1.0 - rhor/ρ_mr)

    M = vs / c_r

    return M - Mach
end


function MachSolver_HelperFunction(Pl::Float64, Pr::Float64,
                                   rhor::Float64, rhol::Float64,
                                   γ_th::Float64, γ_cr::Float64, γ_exp::Float64,
                                   M::Float64, acc_function, etaB::Float64)

    vs = 0.0

    cr = sqrt.( γ_th * Pr / rhor)
    cl = sqrt.( γ_th * Pl / rhol)

    ξ = etaB * acc_function(M)/( 1.0 - acc_function(M))

    if ξ > 0.0
        # solve Density
        f_z(xs) = find_xs(xs, rhor, Pr, Pl, cr, cl,
                            γ_th, γ_cr, γ_exp, ξ)

        xs = find_zero( f_z, 3.0 )

        rho4 = xs*rhor

        # solve pressure
        P34_tot = P4_f(xs, Pr, γ_th, γ_cr, ξ)

        # solve velocity
        v34 = 2.0 * cl / ( γ_th - 1.0 ) * ( 1.0 - (P34_tot / Pl)^γ_exp )

        vs = v34 / ( 1.0 - rhor / rho4 )
    end
    Mach = vs/cr

    return M - Mach
end

function solveMachfromPrCR(Pl::Float64, Pr::Float64,
                     rhor::Float64, rhol::Float64,
                     γ_th::Float64, γ_cr::Float64, γ_exp::Float64,
                     M::Float64, acc_function)

   f(M) = MachSolver_HelperFunction(Pl, Pr,
                                    rhor, rhol,
                                    γ_th, γ_cr, γ_exp,
                                    M, acc_function)

   Mach = find_zero(f, 5.0)

   return Mach
end

function solvePrfromMachCR(Pl::Float64, Pr::Float64,
                     rhor::Float64, rhol::Float64,
                     γ_th::Float64, γ_cr::Float64, γ_exp::Float64,
                     M::Float64, acc_function)

    f(Pr) = MachSolver_HelperFunction(Pl, Pr,
                                     rhor, rhol,
                                     γ_th, γ_cr, γ_exp,
                                     M, acc_function)

    Pr = find_zero(f, 0.03)


    return Pr
end

"""
    Shared functions
"""
function get_ξ(acc_function, M::Float64)
    return ξ = acc_function(M)/(1.0 - acc_function(M))
end

function xs_f(rho4::Float64, rhor::Float64)
    return rho4/rhor
end

function ys_f(xs::Float64, γ_th::Float64, γ_cr::Float64, ξ::Float64)
    return ( ξ * ( xs*(γ_cr - 1.0) - (γ_cr + 1.0) )*xs^γ_th - xs*(γ_th + 1.0) + (γ_th - 1.0) ) /
           ( ξ * ( xs*(γ_cr - 1.0) - (γ_cr + 1.0) )         + xs*(γ_th - 1.0) - (γ_th + 1.0) )
end

function A_f(xs::Float64, γ_th::Float64, γ_cr::Float64, ξ::Float64)
    #return (rho4 - rhor)/(rho4 + rhor)
    return ( ys_f(xs, γ_th, γ_cr, ξ)        + ξ * ( ys_f(xs, γ_th, γ_cr, ξ) - xs^γ_th)         - 1.0 ) /
           ( ys_f(xs, γ_th, γ_cr, ξ) * γ_th + ξ * ( ys_f(xs, γ_th, γ_cr, ξ) - xs^γ_th ) * γ_cr + γ_th )
end

function P_inj_f(xs::Float64, rhor::Float64, Pr::Float64, γ_th::Float64, γ_cr::Float64, ξ::Float64)
    return (( γ_cr - 1.0)/(γ_th - 1.0)) * ξ *
           ( ys_f(xs, γ_th, γ_cr, ξ) - xs^γ_th ) * Pr
end

function P4_f(xs::Float64, Pr::Float64, γ_th::Float64, γ_cr::Float64, ξ::Float64)
    ys(xs) = ys_f(xs, γ_th, γ_cr, ξ)

    return ( ys(xs) + ( γ_cr - 1.0)/(γ_th - 1.0) *
                ξ * ( ys(xs) - xs^γ_th ) )  * Pr
end

function find_xs(xs::Float64, rhor::Float64,
                 Pr::Float64,   Pl::Float64,
                 cr::Float64,   cl::Float64,
                 γ_th::Float64, γ_cr::Float64, γ_exp::Float64,
                 ξ::Float64)

    P4(xs) = P4_f(xs, Pr, γ_th, γ_cr, ξ)
    A(xs) = A_f(xs, γ_th, γ_cr, rhor)

    return ( P4(xs)/Pr - 1.0 ) * A(xs)/( 1.0 + A(xs))  -
              2.0 * γ_th / ( γ_th - 1.0 )^2 *
              cl^2 / cr^2 *
           ( 1.0 - (P4(xs)/Pl)^γ_exp )^2
end

"""
    Functions to solve individual segments of the shocktube
"""


"""
    Pressure
"""
function solveP34(par::SodCRParameters_noCRs, sol::SodCRSolution)
    # Solves the pressure of the middle state in a shocktube
    return P4_f(sol.xs, par.Pr, par.γ_th, par.γ_cr, par.ξ)
end

function solveP2(x::Float64, par::SodCRParameters_noCRs)
    # solves the pressure along the rarefaction wave
    return par.Pl * ( -par.η2 * x / ( par.cl * par.t ) + ( 1.0 - par.η2 ) )^( 1.0/par.γ_exp )
end

function solveP(x::Float64, par::SodCRParameters_noCRs, sol::SodCRSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return par.Pl, # P_tot
               par.Pl, # P_th
               0.0,    # P_cr
               0.0,    # P_cr_p
               0.0     # P_cr_e
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        P2 = solveP2(x, par)
        return P2,
               P2,
               0.0,
               0.0,
               0.0
    elseif -sol.vt * par.t < x <= sol.v34 * par.t
        return sol.P34_tot,
               sol.P34_tot,
               0.0,
               0.0,
               0.0
    elseif sol.v34 * par.t < x <= sol.vs * par.t
        return sol.P34_tot,
               sol.P4_th,
               sol.P4_cr,
               sol.P4_cr * ( 1.0 - par.Pe_ratio),
               sol.P4_cr * par.Pe_ratio
    elseif sol.vs * par.t < x
        return par.Pr,
               par.Pr,
               0.0,
               0.0,
               0.0
    else
        println("Error!")
        return 0.0, 0.0, 0.0, 0.0, 0.0
    end
end


"""
    Density
"""
function solveRho2(x::Float64, par::SodCRParameters_noCRs)
    # solves density along the rarefaction wave
    return par.rhol * ( -par.η2 * x / ( par.cl * par.t ) +
                      ( 1.0 - par.η2 ) )^( 2.0 / (par.γ_th - 1.0 ) )
end

function solveRho3(par::SodCRParameters_noCRs, sol::SodCRSolution)
    # solves density left of contact discontiuity and right of refraction wave
    return par.rhol * ( sol.P34_tot/par.Pl )^( 1.0/par.γ_th )
end


function solveRho4(par::SodCRParameters_noCRs)
    # solves density left of shockfront following
    # Eq. B7-B11 , Pfrommer et. al. 2016

    f_z(xs) = find_xs(xs, par.rhor, par.Pr, par.Pl, par.cr, par.cl,
                        par.γ_th, par.γ_cr, par.γ_exp, par.ξ)

    xs = find_zero( f_z, par.first_guess )

    rho4 = xs*par.rhor

    return rho4

end


function solveRho(x::Float64, par::SodCRParameters_noCRs, sol::SodCRSolution)
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
function solveV34(par::SodCRParameters_noCRs, sol::SodCRSolution)
    # solves velocity along the isopressure region 2-3
    return 2.0 * par.cl / ( par.γ_th - 1.0 ) * ( 1.0 - (sol.P34_tot / par.Pl)^par.γ_exp )
end

function solveV2(x::Float64, par::SodCRParameters_noCRs)
    # solves the velocity along the rarefaction wave
    return ( 1.0 - par.η2 ) * ( x / par.t + par.cl )
end

function solveVs(par::SodCRParameters_noCRs, sol::SodCRSolution)
    # solves the shock velocity
    return sol.v34 / ( 1.0 - par.rhor / sol.rho4 )
end

function solveVt(par::SodCRParameters_noCRs, sol::SodCRSolution)
    # solves the velocity of the tail of the rarefaction wave
    return par.cl - sol.v34/( 1.0 - par.η2 )
end

function solveV(x::Float64, par::SodCRParameters_noCRs, sol::SodCRSolution)
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

function solveSodShockCR_noPrepopulation(x::Array{Float64,1}; par::SodCRParameters_noCRs)

    # set up datatype to store riemann solution
    sol = SodCRSolution(x)

    # transform into rest-frame of contact discontiuity
    x_in = sol.x .- par.x_contact

    # solve Density
    sol.rho4 = solveRho4(par)
    sol.xs  = sol.rho4/par.rhor

    sol.P34_tot = P4_f(sol.xs, par.Pr, par.γ_th, par.γ_cr, par.ξ)

    sol.rho3 = solveRho3(par, sol)

    sol.xr  = sol.rho3/par.rhol

    # solve pressure
    sol.P4_cr = P_inj_f(sol.xs, par.rhor, par.Pr, par.γ_th, par.γ_cr, par.ξ)

    sol.P4_th = sol.P34_tot - sol.P4_cr
    sol.P3_th = sol.P34_tot

    # solve velocity
    sol.v34 = solveV34(par, sol)

    sol.vs = solveVs(par, sol)
    sol.vt = solveVt(par, sol)

    for i=1:length(sol.x)
        sol.v[i] = solveV(x_in[i], par, sol)
    end

    for i = 1:length(sol.x)
        sol.P_tot[i], sol.P_th[i], sol.P_cr[i], sol.P_cr_p[i], sol.P_cr_e[i] = solveP(x_in[i], par, sol)
    end

    for i = 1:length(sol.x)
        sol.rho[i] = solveRho(x_in[i], par, sol)
    end

    sol.Mach = sol.vs/par.cr

    sol.U_tot = sol.P_tot ./ ( (par.γ_th - 1.0) .* sol.rho )
    sol.U_th = sol.P_th ./ ( (par.γ_th - 1.0) .* sol.rho )
    sol.E_cr_p = sol.P_cr_p ./ ( (par.γ_cr - 1.0) .* sol.rho )
    sol.E_cr_e = sol.P_cr_e ./ ( (par.γ_cr - 1.0) .* sol.rho )

    return sol
end

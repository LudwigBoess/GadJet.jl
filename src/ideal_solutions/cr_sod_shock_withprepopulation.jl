using Roots
using Distributions
using NLsolve
using QuadGK
#using ForwardDiff
#using SpecialFunctions
#include("cr_sod_shock_main.jl")
#include("beta_inc_functions.jl")

mutable struct SodCRParameters_withCRs

    rhol::Float64
    rhor::Float64
    Pl::Float64
    Pr::Float64
    Ul::Float64
    Ur::Float64
    P_cr_l::Float64
    P_cr_r::Float64
    E_cr_l::Float64
    E_cr_r::Float64
    cl::Float64
    cr::Float64
    M::Float64
    t::Float64
    x_contact::Float64
    Pe_ratio::Float64
    γ_th::Float64
    γ_cr::Float64
    Δγ::Float64
    γ_exp::Float64
    α::Float64
    β::Float64
    η2::Float64
    # acc_function#::F
    # reacc_function

    function SodCRParameters_withCRs(;rhol::Float64=1.0,  rhor::Float64=0.125,
                                Pl::Float64=0.0,    Pr::Float64=0.0,
                                Ul::Float64=0.0,    Ur::Float64=0.0,
                                P_cr_l::Float64=0.0,      P_cr_r::Float64=0.0,
                                E_cr_l::Float64=0.0,      E_cr_r::Float64=0.0,
                                Mach::Float64=0.0,  t::Float64,
                                x_contact::Float64=70.0,
                                Pe_ratio::Float64=0.01,
                                γ_th::Float64=5.0/3.0,
                                γ_cr::Float64=4.0/3.0,
                                eff_model::Int64=-1)

        γ_exp    = ( γ_th - 1.0 )/( 2.0 * γ_th )
        η2       = (γ_th-1.0)/(γ_th+1.0)
        Δγ       = γ_th - γ_cr

        # calculate Ul and Pl depending on input
        if (Pl == 0.0) & (Ul != 0.0)
            Pl = ( γ_th - 1.0 ) * rhol * Ul
        elseif (Ul == 0.0) & (Pl != 0.0)
            Ul = Pl / ( (γ_th - 1.0) * rhol )
        else
            error("Both Ul and Pl are zero!")
        end

        # calculate B angle dependent efficiency following Pais+ 2018, MNRAS, 478, 5278
        # delta_theta = π/18.0
        # thetaB *= (π/180)
        # etaB = 0.5*( tanh( (theta_crit - thetaB)/delta_theta ) + 1.0 )

        if eff_model == -1
            acc_function = null_acc
        elseif eff_model == 0
            acc_function = KR07
        elseif eff_model == 1
            acc_function = KR13
        elseif eff_model == 2
            acc_function = R_19
        elseif eff_model == 3
            acc_function = CS15
        else
            error("Invalid DSA model selection!\n
                   Pick one of the available models, or solve a pure Hydro shock with:\n
                   SodParameters")

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
            Pr = solvePrfromMach(rhol, rhor, Pl, Mach, γ_th, γ_cr, eff_function)
            Ur = Pr / ( (γ_th - 1.0) * rhor )
        end
        #
        # if Mach == 0.0
        #     Mach = solveMach(Pl, Pr, rhol, rhor, γ_th)
        # end


        cl = sqrt.(γ_eff(Pl, P_cr_l, γ_th, γ_cr) * (Pl + P_cr_l) / rhol)
        cr = sqrt.(γ_eff(Pr, P_cr_r, γ_th, γ_cr) * (Pr + P_cr_r) / rhor)

        α = (γ_cr - 1.0)/(2.0 * Δγ)
        β = (1.0 - γ_th)/(2.0 * Δγ)

        E_cr_r = P_cr_r/( ( γ_cr - 1.0 )*rhor )

        new(rhol, rhor,
            Pl, Pr,
            P_cr_l, P_cr_r,
            Ul, Ur,
            E_cr_r,
            cl, cr,
            Mach, t,
            x_contact,
            Pe_ratio,
            γ_th, γ_cr,
            Δγ,
            γ_exp,
            α, β,
            η2)#,
            #acc_function,
            #acc_function)
    end
end



"""
        Integral function
"""

@inline A(P::Float64, ρ::Float64, γ::Float64) = γ * P * ρ^(-γ)
@inline a(P::Float64, γ::Float64) = γ * P

@inline function x(ρ::Float64, P_th::Float64, P_cr::Float64, par::SodCRParameters_withCRs)
    return  a(P_th, par.γ_th)/(a(P_th, par.γ_th) + a(P_cr, par.γ_cr))
end



@inline function reduced_Beta(x, α, β)
    return 1.0/beta(α, β) * x^(α - 1.0) * ( 1.0 - x )^(β - 1.0)
end

@inline function γ_eff(P_th::Float64, P_cr::Float64, γ_th::Float64, γ_cr::Float64)
    return (γ_cr*P_cr + γ_th*P_th)/(P_th + P_cr)
end

@inline function incomplete_beta(a::Float64, b::Float64, t::Float64)
    return t^(a-1.0) * (1.0-t)^(b-1.0)
end

@inline function I(rho::Float64, P_th::Float64, P_cr::Float64, par::SodCRParameters_withCRs)

    x_ = x(rho, P_th, P_cr,  par)
    A_cr = A(P_cr, rho, par.γ_cr)
    A_th = A(P_th, rho, par.γ_th)
    B(x) = incomplete_beta(par.α, par.β, x)

    result_B = quadgk(B, 0, x_, rtol=1.e-4)

    return √(A_cr)/par.Δγ * (A_cr/A_th)^par.α * result_B
end

"""
    Pressure
"""
function solveP3(par::SodCRParameters_withCRs, sol::SodCRSolution)
    # Solves the pressure of the middle state in a shocktube

    Pcr3 = par.P_cr_l * ( sol.rho3/par.rhol )^par.γ_cr
    P3 = par.Pl * ( sol.rho3/par.rhol )^par.γ_th

    return P3, Pcr3
end

function solveP2(x::Float64; par::SodCRParameters_withCRs)
    # solves the pressure along the rarefaction wave
    return par.Pl * ( -par.η2 * x / ( par.cl * par.t ) + ( 1.0 - par.η2 ) )^( 1.0/par.γ_exp )
end

function solveP(x::Float64, par::SodCRParameters_withCRs, sol::SodCRSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return (par.Pl + par.P_cr_l),
                par.Pl, par.P_cr_l,
                (1.0 - par.Pe_ratio)*par.P_cr_l,
                par.Pe_ratio*par.P_cr_l

    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return 0.0, 0.0, 0.0, 0.0, 0.0

    elseif -sol.vt * par.t < x <= sol.v34 * par.t
        return (sol.P3_th + sol.P3_cr),
                sol.P3_th, sol.P3_cr,
                (1.0 - par.Pe_ratio)*sol.P3_cr,
                par.Pe_ratio*sol.P3_cr

    elseif sol.v34 * par.t < x <= sol.vs * par.t
        return (sol.P4_th + sol.P4_cr),
                sol.P4_th, sol.P4_cr,
                (1.0 - par.Pe_ratio)*par.P4_cr,
                par.Pe_ratio*par.P4_cr

    elseif sol.vs * par.t < x
        return (par.Pr + par.P_cr_r),
                par.Pr, par.P_cr_r,
                (1.0 - par.Pe_ratio)*par.P_cr_r,
                par.Pe_ratio*par.P_cr_r
    else
        println("Error! x = $x  t = $(par.t)")
        return 0.0, 0.0, 0.0, 0.0, 0.0
    end
end


"""
    Density
"""
function solveRho2(x::Float64, par::SodCRParameters_withCRs)
    # solves density along the rarefaction wave

    #f(rho) = I(rho) - I(par.rhol) + x/par.t + sqrt.( A())
    return 0.0

end


function rho34_solver!(F, x, par::SodCRParameters_withCRs)

    xs(rho4)        = rho4/par.rhor
    xr(rho3)        = rho3/par.rhol

    P34(rho3)       = par.P_cr_l*xr(rho3)^par.γ_cr + par.Pl*xr(rho3)^par.γ_th

    Pcr4(rho4)      = par.P_cr_r*xs(rho4)^par.γ_cr

    ε2(rho4, rho3)  = par.E_cr_r*xs(rho4)^par.γ_cr + 1.0/(par.γ_th - 1.0) *
                     ( P34(rho3) - Pcr4(rho4) )

    Pth3(rho3)      = par.Pl * xr(rho3)^par.γ_th
    Pcr3(rho3)      = par.P_cr_l * xs(rho3)^par.γ_cr

    F[1] = (P34(x[1]) - (par.Pr + par.P_cr_r)) * ( xs(x[2]) - 1.0 ) -
            par.rhor * xs(x[2]) * ( I(par.rhol, par.Pl, par.P_cr_l, par) - I(x[1], Pth3(x[1]), Pcr3(x[1]), par) )^2

    F[2] = (P34(x[1]) + (par.Pr + par.P_cr_r)) * ( xs(x[2]) - 1.0 ) +
            2.0 * ( xs(x[2]) * (par.Ur + par.E_cr_r) - ε2(x[2], x[1]) )

end

function solveRho34(par::SodCRParameters_withCRs)

    rho34_helper!(F, x) = rho34_solver!(F, x, par)

    initial_x = [par.rhol, par.rhor]
    nlsolve(rho34_helper!, initial_x)

        #   rho3           rho4
    return initial_x[1], initial_x[2]

end

function solveRho(x::Float64, par::SodCRParameters_withCRs, sol::SodCRSolution)
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
function solveV34(par::SodCRParameters_withCRs, sol::SodCRSolution)
    # solves velocity along the isopressure region 2-3
    return sqrt.( (sol.P34_tot - (par.Pr + par.P_cr_r ) ) * (sol.rho4 - par.rhor)/sol.rho4*par.rhor )
end

function solveV2(x::Float64, par::SodCRParameters_withCRs)
    # solves the velocity along the rarefaction wave
    return 0.0
end

function solveVs(par::SodCRParameters_withCRs, sol::SodCRSolution)
    # solves the shock velocity
    return sol.v34 * sol.rho4 / ( sol.rho4 - par.rhor )
end

function solveVt(par::SodCRParameters_withCRs, sol::SodCRSolution)
    # solves the velocity of the tail of the rarefaction wave
    return I(sol.rho3, sol.P3_th, sol.P3_cr, par) - I(par.rhol, par.Pl, par.P_cr_l, par) +
            sqrt.(
                A(sol.P3_cr, sol.rho3, par.γ_cr)*sol.rho3^(par.γ_cr - 1.0) +
                A(sol.P3_th, sol.rho3, par.γ_th)*sol.rho3^(par.γ_th - 1.0) )
end

function solveV(x::Float64, par::SodCRParameters_withCRs, sol::SodCRSolution)
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


function solveSodShockCR_withPrepopulation(x::Array{Float64,1}; par::SodCRParameters_withCRs)

    # set up datatype to store riemann solution
    sol = SodCRSolution(x)

    # transform into rest-frame of contact discontiuity
    x_in = sol.x .- par.x_contact

    # solve Pressure
    sol.rho3, sol.rho4 = solveRho34(par)


    sol.xs  = sol.rho4/par.rhor
    sol.xr  = sol.rho3/par.rhol

    sol.P34_tot = par.P_cr_l*sol.xr^par.γ_cr + par.Pl*sol.xr^par.γ_th

    sol.P3_cr = par.P_cr_l*(sol.rho3/par.rhol)^par.γ_cr
    sol.P3_th = sol.P34_tot - sol.P3_cr

    sol.P4_cr = par.P_cr_r*sol.xs^par.γ_cr
    sol.P4_th = sol.P34_tot - sol.P4_cr

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

# x_in = collect(50.0:0.01:100.0)
#
# Pl = 63.499
# Pr = 0.1
# Pcrl = 0.3*Pl
# Pcrr = 0.3*Pr
# par = SodCRParameters_withCRs(Pl=Pl, Pr=Pr, P_cr_l=Pcrl, P_cr_r=Pcrr, t=1.5, eff_model=-1)
#
# sol = solveSodShockCR_withPrepopulation(x_in, par=par)
#
# B(x) = incomplete_beta(par.α, par.β, x)
#
# result_B = quadgk(B, 0, 4.74, rtol=1.e-4)

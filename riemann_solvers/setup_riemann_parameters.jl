include("sod_shock.jl")
include("cr_sod_shock_main.jl")
include("cr_sod_shock_noprepopulation.jl")
include("cr_sod_shock_withprepopulation.jl")

function RiemannParameters(;rhol::Float64=1.0,        rhor::Float64=0.125,
                            Pl::Float64=0.0,          Pr::Float64=0.0,
                            Ul::Float64=0.0,          Ur::Float64=0.0,
                            P_cr_l::Float64=0.0,      P_cr_r::Float64=0.0,
                            E_cr_l::Float64=0.0,      E_cr_r::Float64=0.0,
                            Bl::Array{Float64,1} = zeros(3),
                            Br::Array{Float64,1} = zeros(3),
                            Mach::Float64=0.0,        t::Float64,
                            x_contact::Float64=70.0,
                            γ_th::Float64=5.0/3.0,
                            Pe_ratio::Float64=0.01,
                            γ_cr::Float64=4.0/3.0,
                            thetaB::Float64=0.0,
                            theta_crit::Float64=(π/4.0),
                            dsa_model::Int64=-1)


    # Error handling
    if (Pl == 0.0 && Ul == 0.0)
        error("Both Ul and Pl are zero!")
    end # error handling


    # pure hydro case -> no pre-existing CRs and no DSA model selected
    if ( (P_cr_l == 0.0 && P_cr_r == 0.0) ||
         (E_cr_l == 0.0 && E_cr_r == 0.0) ) &&
         dsa_model == -1

         println("Setting up parameters for pure hydro Sod-shock.")

         return SodParameters(rhol, rhor, Pl, Pr,
                              Ul, Ur, Mach, t,
                              x_contact, γ_th)

    elseif dsa_model != -1

        println("Setting up parameters for Sod-shock with CR acceleration.")

        if ( (P_cr_l == 0.0 && P_cr_r == 0.0) ||
             (E_cr_l == 0.0 && E_cr_r == 0.0) )

            println("No seed CRs.")
            return SodCRParameters_noCRs(rhol=rhol, rhor=rhor,
                                         Pl=Pl, Pr=Pr, Ul=Ul, Ur=Ur,
                                         Mach=Mach, t=t, x_contact=x_contact,
                                         Pe_ratio=Pe_ratio, γ_th=γ_th, γ_cr=γ_cr,
                                         thetaB=thetaB, theta_crit=theta_crit,
                                         dsa_model=dsa_model)
        else
            println("With seed CRs.")
            return SodCRParameters_withCRs(rhol=rhol, rhor=rhor,
                                           Pl=Pl, Pr=Pr, Ul=Ul, Ur=Ur,
                                           P_cr_l=P_cr_l, P_cr_r=P_cr_r,
                                           E_cr_l=E_cr_l, E_cr_r=E_cr_r,
                                           Mach=Mach, t=t, x_contact=x_contact,
                                           Pe_ratio=Pe_ratio, γ_th=γ_th, γ_cr=γ_cr,
                                           thetaB=thetaB, theta_crit=theta_crit,
                                           dsa_model=dsa_model)
        end #
    end # dsa_model != -1

    if ( (P_cr_l != 0.0 && P_cr_r != 0.0) ||
         (E_cr_l != 0.0 && E_cr_r != 0.0) ) &&
         dsa_model == -1

         println("Setting up parameters for multicomponent shock without CR acceleration.")
         return SodCRParameters_withCRs(rhol, rhor, Pl, Pr,
                                         Ul, Ur, P_cr_l, P_cr_r,
                                         E_cr_l, E_cr_r,
                                         Mach, t, x_contact,
                                         Pe_ratio, γ_th, γ_cr,
                                         eff_model)

    end # Seed CRs without acc

end # RiemannParameters


"""
    Multiple dispatch for solve function
"""
# Pure hydro Sod-shock
solve(x::Array{Float64,1}, par::SodParameters) = solveSodShock(x, par=par)
solve(x::Array{Float32,1}, par::SodParameters) = solveSodShock(Float64.(x), par=par)

# CR Sod shock
solve(x::Array{Float64,1}, par::SodCRParameters_noCRs)   = solveSodShockCR_noPrepopulation(x, par=par)
solve(x::Array{Float64,1}, par::SodCRParameters_noCRs)   = solveSodShockCR_noPrepopulation(Float64.(x), par=par)

solve(x::Array{Float64,1}, par::SodCRParameters_withCRs) = solveSodShockCR_withPrepopulation(x, par=par)
solve(x::Array{Float32,1}, par::SodCRParameters_withCRs) = solveSodShockCR_withPrepopulation(Float64.(x), par=par)

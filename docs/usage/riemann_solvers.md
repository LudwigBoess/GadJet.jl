Riemann Solvers
===============

GadJet.jl provides a number of exact riemann solvers.
So far these are for
  * Sod shock, pure hydro
  * Sod shock, with CR acceleration

Setup
-----

To get the exact solution to a Sod shock you first need to set up the inital conditions.
You can do this with the helper function `RiemannParameters` that contains all parameters for all possible configurations:

```julia
RiemannParameters(;rhol::Float64=1.0, rhor::Float64=0.125,      # density left and right (L&R)
                   Pl::Float64=0.0,   Pr::Float64=0.0,          # pressure L&R
                   Ul::Float64=0.0,   Ur::Float64=0.0,          # internal energy L&R
                   P_cr_l::Float64=0.0, P_cr_r::Float64=0.0,    # CR pressure L&R
                   E_cr_l::Float64=0.0, E_cr_r::Float64=0.0,    # CR energy L&R
                   Bl::Array{Float64,1} = zeros(3),             # B-field left
                   Br::Array{Float64,1} = zeros(3),             # B-field right
                   Mach::Float64=0.0,                           # target Mach number
                   t::Float64,                                  # time of the solution
                   x_contact::Float64=70.0,                     # position of the contact discontinuity along the tube
                   γ_th::Float64=5.0/3.0,                       # adiabatic index of the gas               
                   γ_cr::Float64=4.0/3.0,                       # adiabatic index of CRs
                   Pe_ratio::Float64=0.01,                      # ratio of proton to electron energy in acceleration
                   thetaB::Float64=0.0,                         # angle between magnetic field and shock normal
                   theta_crit::Float64=(π/4.0),                 # critical angle for B/Shock angle efficiency
                   dsa_model::Int64=-1,                         # diffuse shock acceleration model
                   xs_first_guess::Float64=4.7)                 # first guess of the resulting shock compression
```

To set up a standard Sod shock you need to supply it with pressure/energy values for left and right state, or with pressure/energy values for the left state and a target Mach number.

A minimal working version would be, for a shock with Mach 10, at time = 1.5:

```julia
par = RiemannParameters(Ul=100.0, Mach=10.0, t=1.5)
```

This returns a parameter object for a pure hydro Sod shock:

```julia
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

end
```

A minimal working version for the solution of the CR shock discussed in Pfrommer+16 (doi:10.1093/mnras/stw2941) would be:

```julia
par = RiemannParameters(Pl=63.499, Pr=0.1, t=1.5, dsa_model=4)
```

This also returns a parameter object: `SodCRParameters_noCRs` which can be found in cr_sod_shock_noprepopulation.jl .

Solving the shock
-----------------

To solve the shock with the given initial condition you just need to call

```julia
sol = solve(x, par)
```

with par being either of the above mentioned parameter objects, multiple dispatch will take care of the rest.

`x` has to be an array with either sample positions along the tube, or your actual particle positions, to make calculating errors easier. You can also just pass it an array with a single position, if you're only interested in that specific part of the shock ( e.g. `x = [86.0]` for the center of the postshock region.)

This will return a solution object depending on which shock you're solving.

For the pure hydro case this is:

```julia
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
end
```

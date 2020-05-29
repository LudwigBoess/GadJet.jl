using ProgressMeter
using Interpolations
using Roots

struct Hernquist

    Npart::Int64
    Mtot::Float64
    r0::Float64
    rmax::Float64
    rho0::Float64

    function Hernquist(Npart::Int64, Mtot::Float64, r0::Float64, rmax::Float64)

        mass_analytic(r, ρ0, r0) = 2π * ρ0 * r0^3 * (r/r0)^2 / ( 1.0 + r/r0 )^2
        root_helper(ρ0) = mass_analytic(rmax, ρ0, r0) - Mtot

        ρ0_ideal = find_zero(root_helper, 1.0)

        new(Npart, Mtot, r0, rmax, rho0)

    end
end

## Analytic profiles

function density_profile(r::Float64, halo::Hernquist)
    x = r/halo.r0
    return halo.rho0/( x * ( 1.0 + x)^3 )
end

function culm_mass_profile(r::Float64, halo::Hernquist)
    x = r/halo.r0
    return 2π * halo.rho0 * halo.r0^3 * x^2 / ( 1.0 + x )^2
end

# functions to sample positions
function sample_xr_xv(halo::Hernquist)

    xr = 0.0
    xv = 0.0
    e = 10.0
    while ( e > 0.0 || e < (-halo.Mtot/halo.r0) )

        # sample r from system size
        xr = halo.rmax * rand()
        # sample v from escape velocity
        xv = sqrt(2.0*halo.Mtot/halo.r0) * rand()

        # Make sure energy is negative -> particle is bound!
        e = 0.5xv^2 - halo.Mtot/(xr+halo.r0)
    end

    return xr, xv, e
end


function get_pn(halo::Hernquist, xr::Float64, xv::Float64, e::Float64)

    q = sqrt(-1.0*halo.r0*e/halo.Mtot)
    q2 = q*q

    fq = ( 3.0*asin(q) + q*sqrt(1.0 - q2) *
         ( 1.0 - 2.0 * q2) * (8.0*q2*q2 - 8.0*q2 - 3.0)) /
         ( 1.0 - q2 )^2.5

    return ( xr*xr/(halo.r0*halo.r0) ) * ( xv*xv / ( halo.Mtot/halo.r0 )) * fq / 1.20223157581242
end

function sample_particle(xr::Float64, xv::Float64)

    # sample position
    cth   = 2.0 * ( rand() - 0.5 )
    signs = 2.0 * ( rand() - 0.5 )
    sth   = signs/abs(signs) * sqrt( 1.0 - cth*cth)
    phi   = 2π * rand()

    x = xr*sth*cos(phi)
    y = xr*sth*sin(phi)
    z = xr*cth

    # sample velocity
    cth   = 2.0 * ( rand() - 0.5 )
    signs = 2.0 * ( rand() - 0.5 )
    sth   = signs/abs(signs) * sqrt( 1.0 - cth*cth)
    phi   = 2π * rand()

    vx = xv*sth*cos(phi)
    vy = xv*sth*sin(phi)
    vz = xv*cth

    return x, y, z, vx, vy, vz
end


function sample_particle_pos_vel(halo::Hernquist)

    x  = zeros(halo.Npart)
    y  = zeros(halo.Npart)
    z  = zeros(halo.Npart)

    vx = zeros(halo.Npart)
    vy = zeros(halo.Npart)
    vz = zeros(halo.Npart)

    P = Progress(halo.Npart)
    idx = 0

    Np = 1
    while Np < halo.Npart

        xr, xv, e = sample_xr_xv(halo)

        if (e > 0.0)
            error("e = $e")
        end

        pn = get_pn(xr, xv, e, halo)

        if (rand() <= pn)
            x[Np], y[Np], z[Np], vx[Np], vy[Np], vz[Np] = sample_particle(xr, xv)

            Np  += 1
            idx += 1
            ProgressMeter.update!(P, idx)

        end
    end

    return [x y z], [vx vy vz]
end

function calculate_particle_masses(halo::Hernquist)
    mk = halo.Mtot * halo.rmax * halo.rmax / (halo.rmax + halo.r0)^2 / halo.Npart
    return mk .* ones(halo.Npart)
end
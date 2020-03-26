using Unitful
using UnitfulAstro

# set up proton number density unit
@unit n_p "N_p/cm^3" ProtonNumberDensity 1u"mp/cm^3" true
@unit n_e "N_e/cm^3" ElectronNumberDensity 1u"me/cm^3" true

# needed by unitful
const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(n_p)
    Unitful.register(n_e)
end

struct GadgetPhysicalUnits

    x_cgs::typeof(1u"cm")         # position in cm
    v_cgs::typeof(1u"cm/s")       # velocity in cm/s
    m_cgs::typeof(1u"g")          # mass in g

    t_s::typeof(1u"s")            # time in sec
    t_Myr::typeof(1u"Myr")        # time in Myr

    E_cgs::typeof(1u"cgs")        # energy in erg
    E_eV::typeof(1u"eV")          # energy in eV

    B_cgs::typeof(1u"Gs")         # magnetic field in Gauss

    rho_cgs::typeof(1u"g/cm^3")        # density in g/cm^3
    rho_ncm3::typeof(1u"N_p/cm^3")     # density in N_p/cm^3

    T_K::typeof(1u"K")                 # temperature in K

    P_th_cgs::typeof(1u"Ba")           # pressure in
    #P_th_si::typeof(1u"Bar")
    P_CR_cgs::typeof(1u"Ba")
    #P_CR_si::typeof(1u"Bar")

    function GadgetPhysicalUnits(l_unit=3.085678e21u"cm", m_unit=1.989e43u"g", v_unit=1.e5u"cm/s";
                                 a_scale=1.0, hpar=1.0, γ_th=5.0/3.0, γ_CR=4.0/3.0, xH=0.76)

        x_cgs  = l_unit * a_scale/hpar
        v_cgs  = v_unit * sqrt(a_scale)
        m_cgs  = m_unit * hpar
        t_unit = l_unit / v_unit
        t_s    = t_unit * sqrt(a_scale) / hpar   # in sec
        t_Myr  = t_s |> u"Myr"



        E_cgs = m_cgs * v_cgs^2 |> u"erg"
        E_eV = E_cgs |> u"eV"

        # B_unit = 1.0    # gadget outputs in cgs
        B_cgs = 1.0u"Gs"    # gadget outputs in cgs
        # B_si = 1.e-4

        rho_cgs = m_unit/l_unit^3 * hpar^2 / a_scale^3
        rho_ncm3 = rho_cgs |> u"n_p"

        yhelium = ( 1.0 - xH ) / ( 4.0 * xH )
        mean_mol_weight = (1.0 + 4.0 * yhelium) / (1.0 + 3.0 * yhelium + 1.0)

        T_cgs = (γ_th - 1.0) * v_unit^2 * 1u"mp" * mean_mol_weight / 1u"k" |> u"K"

        P_th_cgs = (1 + z)^(3*γ_th) * E_unit / l_unit^3 * hpar^2  |> u"Ba"
        P_CR_cgs = (1 + z)^(3*γ_CR) * E_unit / l_unit^3 * hpar^2  |> u"Ba"

        new(x_cgs, v_cgs, m_cgs,
            t_s, t_Myr,
            E_cgs, E_eV,
            B_cgs,
            rho_cgs, rho_ncm3,
            T_cgs,
            P_th_cgs,
            P_CR_cgs,
            )

    end
end

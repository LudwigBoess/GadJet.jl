struct GadgetPhysicalUnits

    x_cgs::Float64          # position in cm
    v_cgs::Float64          # velocity in cm/s
    m_cgs::Float64          # mass in g

    t_s::Float64            # time in sec
    t_Myr::Float64          # time in Myr

    E_cgs::Float64          # energy in erg
    E_eV::Float64           # energy in eV

    # B_unit::Float64
    # B_cgs::Float64

    rho_cgs::Float64        # density in g/cm^3
    rho_ncm3::Float64       # density in N/cm^3

    T_K::Float64            # temperature in K
    #
    # P_th_cgs::Float64
    # P_th_si::Float64
    # P_CR_cgs::Float64
    # P_CR_si::Float64

    function GadgetPhysicalUnits(l_unit=3.085678e21, m_unit=1.989e43, v_unit=1.e5;
                                 a_scale=1.0, hpar=1.0, γ_th=5.0/3.0, γ_CR=4.0/3.0, xH=0.76)

        x_cgs  = l_unit * a_scale/hpar
        v_cgs  = v_unit * sqrt(a_scale)
        m_cgs  = m_unit * hpar
        t_unit = l_unit / v_unit
        t_s    = t_unit * sqrt(a_scale) / hpar   # in sec
        t_Myr  = t_unit/3.1536e+13   # in Myr



        E_cgs = m_cgs * v_cgs^2
        E_eV = E_cgs * 6.242e+11

        # B_unit = 1.0    # gadget outputs in cgs
        # B_cgs = 1.0    # gadget outputs in cgs
        # B_si = 1.e-4

        # proton mass in g
        mp = 1.6726219e-24

        rho_cgs = m_unit/l_unit^3 * hpar^2 / a_scale^3
        rho_ncm3 = rho_cgs/mp

        yhelium = ( 1.0 - xH ) / ( 4.0 * xH )
        mean_mol_weight = (1.0 + 4.0 * yhelium) / (1.0 + 3.0 * yhelium + 1.0)
        kb = 1.380658e-16  # Boltzmann constant

        T_cgs = (γ_th - 1.0) * v_unit^2 * mp * mean_mol_weight / kb

        # P_th_cgs = (1 + z)^(3*γ_th) * E_unit / l_unit^3 * hpar^2
        # P_CR_cgs = (1 + z)^(3*γ_CR) * E_unit / l_unit^3 * hpar^2

        new(x_cgs, v_cgs, m_cgs,
            t_s, t_Myr,
            E_cgs, E_eV,
#            B_unit, B_cgs, B_si,
            rho_cgs, rho_ncm3,
            T_cgs
            # P_th_cgs,
            # P_CR_cgs,
            )

    end

end

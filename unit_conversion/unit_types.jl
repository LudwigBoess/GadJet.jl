mutable struct GadgetUnitFactors

    l_unit::Float64
    m_unit::Float64
    v_unit::Float64

    z::Float64
    hpar::Float64

    γ_th::Float64
    γ_CR::Float64

    t_unit::Float64
    t_Myr::Float64

    E_unit::Float64
    E_cgs::Float64
    E_si::Float64
    E_eV::Float64

    B_unit::Float64
    B_cgs::Float64
    B_si::Float64

    rho_unit::Float64
    rho_cgs::Float64

    T_cgs::Float64

    P_th_cgs::Float64
    P_th_si::Float64
    P_CR_cgs::Float64
    P_CR_si::Float64

    function GadgetUnitFactors(l_unit=3.085678e21, m_unit=1.989e43, v_unit=1.e5;
                         z=0., hpar=0.7, γ_th=5.0/3.0, γ_CR=4.0/3.0, xH=0.76)

        t_unit = l_unit/v_unit
        t_Myr = t_unit/3.1536e+13   # in Myr

        E_unit = m_unit*v_unit^2
        E_cgs = 1.0/E_unit
        E_si = 1.e-7 * E_cgs
        E_eV = E_cgs * 6.242e+11

        B_unit = 1.0    # gadget outputs in cgs
        B_cgs = 1.0    # gadget outputs in cgs
        B_si = 1.e-4

        rho_unit = m_unit/l_unit^3  # 10 M_⊙/pc^3
        rho_cgs = (1.0 + z)^3 * m_unit / l_unit^3 * hpar^2


        yhelium = ( 1.0 - xH ) / ( 4.0 * xH )
        mean_mol_weight = (1.0 + 4.0 * yhelium) / (1.0 + 3.0 * yhelium + 1.0)
        bk = 1.380658e-16
        prtn = 1.672623e-24

        T_cgs = (γ_th - 1.0) * v_unit^2 * prtn * mean_mol_weight / bk

        P_th_cgs = (1 + z)^(3*γ_th) * E_unit / l_unit^3 * hpar^2
        P_th_si = 0.1*P_th_cgs
        P_CR_cgs = (1 + z)^(3*γ_CR) * E_unit / l_unit^3 * hpar^2
        P_CR_si = 0.1*P_CR_cgs


        new(l_unit, m_unit, v_unit, z, hpar, γ_th, γ_CR,
            t_unit, t_Myr,
            E_unit, E_cgs, E_si, E_eV,
            B_unit, B_cgs, B_si,
            rho_unit, rho_cgs,  # rho_si!
            T_cgs,
            P_th_cgs, P_th_si,
            P_CR_cgs, P_CR_si
            )

    end

end

"""
    This implementation is based on the port from GSL to Python by
    J. Michael Burgess for his pynchrotron package:
    https://github.com/grburgess/pynchrotron
"""

GSL_SQRT_DBL_EPSILON = 1.4901161193847656e-08
GSL_LOG_DBL_MIN = -7.0839641853226408e02

global const synch_c0 = π / sqrt(3.0)

global const c01 = 0.2257913526447274323630976
global const cond1 = 2 * sqrt(2.0) * GSL_SQRT_DBL_EPSILON
global const cond3 = -8.0 * GSL_LOG_DBL_MIN / 7.0

global const synchrotron1_data = [  30.364682982501076273,
                                    17.079395277408394574,
                                    4.560132133545072889,
                                    0.549281246730419979,
                                    0.372976075069301172e-01,
                                    0.161362430201041242e-02,
                                    0.481916772120371e-04,
                                    0.10512425288938e-05,
                                    0.174638504670e-07,
                                    0.22815486544e-09,
                                    0.240443082e-11,
                                    0.2086588e-13,
                                    0.15167e-15
                                  ]

global const synchrotron1_cs_order = 12
global const synchrotron1_cs_a = -1.0
global const synchrotron1_cs_b = 1.0

global const synchrotron1a_data = [ 2.1329305161355000985,
                                    0.741352864954200240e-01,
                                    0.86968099909964198e-02,
                                    0.11703826248775692e-02,
                                    0.1645105798619192e-03,
                                    0.240201021420640e-04,
                                    0.35827756389389e-05,
                                    0.5447747626984e-06,
                                    0.838802856196e-07,
                                    0.13069882684e-07,
                                    0.2053099071e-08,
                                    0.325187537e-09,
                                    0.517914041e-10,
                                    0.83002988e-11,
                                    0.13352728e-11,
                                    0.2159150e-12,
                                    0.349967e-13,
                                    0.56994e-14,
                                    0.9291e-15,
                                    0.152e-15,
                                    0.249e-16,
                                    0.41e-17,
                                    0.7e-18
                                ]

global const synchrotron1a_cs_order = 22
global const synchrotron1a_cs_a = -1.0
global const synchrotron1a_cs_b = 1.0

global const synchrotron2_data = [  0.4490721623532660844,
                                    0.898353677994187218e-01,
                                    0.81044573772151290e-02,
                                    0.4261716991089162e-03,
                                    0.147609631270746e-04,
                                    0.3628633615300e-06,
                                    0.66634807498e-08,
                                    0.949077166e-10,
                                    0.1079125e-11,
                                    0.10022e-13,
                                    0.77e-16,
                                    0.5e-18,
                                ]


global const synchrotron2_cs_order = 11
global const synchrotron2_cs_a = -1.0
global const synchrotron2_cs_b = 1.0



@inline function cheb_eval(coeff::typeof([1.0]), order::Integer, a::Real, b::Real, x::Real)

    d = 0.0
    dd = 0.0

    y = (2x - a - b) / (b - a)
    y2 = 2y

    @inbounds for j = (order+1):-1:2
        temp = d
        d = y2 * d - dd + coeff[j]
        dd = temp
    end

    temp = d
    d = y * d - dd + 0.5 * coeff[1]

    return d
end


@inline function synchrotron_kernel(x::Real)

    if x < cond1

        z = x^(1/3)
        cf = 1 - 8.43812762813205e-01 * z * z
        return 2.14952824153447863671 * z * cf

    elseif x <= 4.0

        px = x^(1/3)
        px11 = px^11
        t = x^2 / 8.0 - 1.0
        result_c1 = cheb_eval( synchrotron1_data,
                               synchrotron1_cs_order,
                               synchrotron1_cs_a,
                               synchrotron1_cs_b,
                               t )

        result_c2 = cheb_eval( synchrotron2_data,
                               synchrotron2_cs_order,
                               synchrotron2_cs_a,
                               synchrotron2_cs_b,
                               t )

        return px * result_c1 - px11 * result_c2 - synch_c0 * x

    elseif x < cond3

        t = (12.0 - x) / (x + 4.0)

        result_c1 = cheb_eval( synchrotron1a_data,
                               synchrotron1a_cs_order,
                               synchrotron1a_cs_a,
                               synchrotron1a_cs_b,
                               t )

        return sqrt(x) * result_c1 * exp(c01 - x)

    else
        return 0.0
    end
end


@inline function ndensity_integral(bound_low::Real, bound_up::Real, norm::Real,
                                  slope::Real, density_in::Real)

    nb = 4π * norm * bound_low^3 / density_in
    density = nb * ( (bound_up/bound_low)^(3.0-slope) - 1.0 ) / ( 3.0 - slope )

	Δq = 1.e-6 # "softening" to avoid divergence in denominator
    if ( (3.0 - Δq) < slope < (3.0 + Δq) )

	    slope_var = ( slope - 3.0 ) / Δq
	    density2 = nb * log(bound_up/bound_low)
	    if ( slope_var != 0.0 )
	      density = density * slope_var + density2 * (1.0 - slope_var)
	    else
	      density = density2
	  	end
	end

	density
end


function calculate_synch_intensity(CReNorm, CReSlope, bounds, bin_width::Real,
                                   B::Real, density::Real, ν0::Real=1.4e9)

    m_e = 9.10953e-28
    c_l = 2.9979e10
    qe  = 4.8032e-10
    erg2eV = 6.242e+11

    E_cntr_prefac = m_e*c_l*c_l * sqrt(2.0/3.0 * (2π * m_e*c_l)/qe /0.29)
    nu_c_prefac = 3.0 * qe / (4π * m_e^3 * c_l^5)
    j_nu_prefac = qe * qe * qe * sqrt(3.0) / (m_e * c_l * c_l)

    LMB_SPECTRAL_CRs = length(CReNorm)

    F     = zeros(LMB_SPECTRAL_CRs)
    F_mid = zeros(LMB_SPECTRAL_CRs)
    E     = zeros(LMB_SPECTRAL_CRs)
    dE    = zeros(LMB_SPECTRAL_CRs)
    J     = zeros(LMB_SPECTRAL_CRs)
    N     = zeros(LMB_SPECTRAL_CRs)
    N_mid = zeros(LMB_SPECTRAL_CRs)
    x     = zeros(LMB_SPECTRAL_CRs)

    E_min = bounds[1] * m_e * c_l^2
    E_max = bounds[end] * m_e * c_l^2
    di = log(E_max/E_min)/LMB_SPECTRAL_CRs


    for i = 1:LMB_SPECTRAL_CRs
        E[i]     = bounds[i] * m_e * c_l^2
        x[i]     = ν0 / (nu_c_prefac * E[i]^2 * B)

		# old version
        # N[i]     = CReNorm[i]
        # bound_up = bounds[i] * 10^(bin_width*0.5)
        # N_mid[i] = CReNorm[i] * (bound_up/bounds[i])^(-CReSlope[i])

		# integrate over half a bin
		bound_up = bounds[i] * 10^(bin_width*0.5)
		N[i]     = ndensity_integral(bounds[i], bound_up, CReNorm[i],
									CReSlope[i], density)

		Norm_mid = CReNorm[i] * (bound_up/bounds[i])^(-CReSlope[i])
        N_mid[i] = ndensity_integral(bound_up, bounds[i+1], Norm_mid,
									CReSlope[i], density)
    end

    dE[1] = E[1] - E_min * exp(-1*di)
    for i = 2:LMB_SPECTRAL_CRs
        dE[i] = E[i] - E[i-1]
    end

    for i = 2:LMB_SPECTRAL_CRs

        K = synchrotron_kernel(x[i])
        F[i] = N[i] * K

        x_mid = 0.5 * (x[i-1] + x[i])
        K_mid = synchrotron_kernel(x_mid)

        F_mid[i] = N_mid[i] * K_mid

        J[i] = dE[i] / 6.0 * (F[i] + F[i-1] + 4*F_mid[i])
    end

    return J
end

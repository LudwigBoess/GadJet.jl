e_unit = 1.989e53 # [erg]
m_unit = 1.989e43 # [g]
l_unit = 3.085678e21 # [cm]
v_unit = 1.e5 # [cm/s]
mp = 1.6726219e-24 # [g]

m_unit*v_unit^2

function gaussToGadget(Gauss; GU)

    g = Gauss * sqrt.(GU.l_unit^3/GU.e_unit)

    return g

end

function secToGadget()

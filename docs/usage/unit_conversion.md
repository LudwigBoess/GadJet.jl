Unit Conversion
===============

!!! This is not finished !!!

You can convert the internal units of Gadget into cgs units by defining the object `GadgetPhysicalUnits`:

```julia
GadgetPhysicalUnits(l_unit=3.085678e21, m_unit=1.989e43, v_unit=1.e5;
                    a_scale=1.0, hpar=1.0,
                    γ_th=5.0/3.0, γ_CR=4.0/3.0, xH=0.76)
```

This returns an object of type `GadgetPhysicalUnits` with the following properties:

```julia
struct GadgetPhysicalUnits

    x_cgs::Float64          # position in cm
    v_cgs::Float64          # velocity in cm/s
    m_cgs::Float64          # mass in g

    t_s::Float64            # time in sec
    t_Myr::Float64          # time in Myr

    E_cgs::Float64          # energy in erg
    E_eV::Float64           # energy in eV

    rho_cgs::Float64        # density in g/cm^3
    rho_ncm3::Float64       # density in N/cm^3

    T_K::Float64            # temperature in K
end
```

To convert, say positions of gas particles from a cosmological simulation to physical units you can use:

```julia

h     = read_header(filename)

pos   = read_snapshot(filename, "POS", parttype=0)

GU    = GadgetPhysicalUnits(a_scale=h.time, hpar=h.h0)

pos .*= GU.x_cgs

```

If you have different units than the standard Gadget ones you can call the object cunstructor with different values

```julia
GU = GadgetPhysicalUnits(your_l_unit, your_m_unit, your_v_unit; kwargs...)
```

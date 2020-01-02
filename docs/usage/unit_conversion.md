Unit Conversion
===============

!!! This is not finished !!!

You can convert the internal units of Gadget into cgs units by defining the object `GadgetUnitFactors`:

```julia
GU = GadgetUnitFactors(l_unit=3.085678e21, m_unit=1.989e43, v_unit=1.e5;
                       z=nothing, cosmo_a=nothing, hpar=0.7,
                       γ_th=5.0/3.0, γ_CR=4.0/3.0, xH=0.76)
```

SPH mapping
===========


Internal Module
---------------


External Programs
-----------------
GadJet.jl provides helper function for two external sph mapping Codes: Smac and P-Smac2.

### P-Smac2
P-Smac2 by Julius Donnert (https://github.com/jdonnert/Smac2) is an advanced mapping code for a multitude of different quantities. To run a mapping and plotting loop from a Julia script you need to update the parameter files on the fly.
The function `write_smac2_par` provides this functionality.

```julia
write_smac2_par(x, y, z,
                euler_angle_0, euler_angle_1, euler_angle_2,
                xy_size, z_depth, xy_pix::Int64,
                input_file, output_file, path,
                effect_module::Int64=0, effect_flag::Int64=0)
```

### Smac
Smac isn't public unfortunately. So these functions are just for my personal use.

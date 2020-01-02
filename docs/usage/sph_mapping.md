SPH mapping
===========


Internal Module
---------------
You can map SPH data to a grid using the function:

```julia
function sphMapping(Pos, HSML, M, ρ, Bin_Quant;
		            param::mappingParameters, kernel,
		            show_progress::Bool=true,
				    # conserve_quantities::Bool=false,  # not used yet
				    dimensions::Int=2)

	[...]

end
```

At the current state this computes the grid pixels by evaluating the kernel at the center of the pixel, like the basic version of SPLASH by Daniel Price.

### Setup
To map the data you need to define the mapping parameters via the `mappingParameters` object:

```julia
par = mappingParameters(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax],
						Npixels=200)
```
Instead of Npixels you can also give the keyword argument `pixelSideLength` if you prefer to define your image that way.

You also need to choose the kernel you used in the simulation. I implemented the following ones:

```julia
k = Cubic()
k = Quintic()
k = WendlandC4()
k = WendlandC6()
```

### Mapping
With the setup done you can now map (e.g.) density of your data using the function above as:

```julia
image = sphMapping(x, hsml, m, rho, rho, param=par, kernel=k)
```

Replacing the second `rho` with any other quantity would map that quantity of course.
Please note: This function doesn't do any unit conversion for you, so you need to convert to the desired units beforehand. See the chapter on unit conversion for usage.

Image now contains a 2D array with the binned data and can easily be plotted with `imshow()` from any plotting package of your choosing.

!!! I'm not very happy with the performance of this module! If you want to contribute to making this faster, your help would be very much appreciated! Please feel free to contact me via the email adress provided in "Contributing" !!!



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

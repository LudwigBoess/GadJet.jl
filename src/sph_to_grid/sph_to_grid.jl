import GR
using Base.Threads

"""
glimpse(filename, blockname[, ... ])

	Gives a glimpse at the simulation by mapping the sph data.
"""
function glimpse(filename::String, blockname::String,
                 center_pos::Array{Float64,1}=[123456.7, 123456.7, 123456.7],
                 dx::Float64=0.0, dy::Float64=0.0, dz::Float64=0.0;
			     kernel_name::String="WC6",
                 resolution::Int64=500, run_dummy::Bool=true,
				 conserve_quantities::Bool=false,
			     verbose::Bool=true, plot::Bool=true)

    if verbose
		@info "Reading data..."
    end
    # read header of snapshot
    h = head_to_obj(filename)
    info = read_info(filename, verbose=false)

    if info == 1
        error("Glimpse only works for snapshots with info block, sorry!")
    end

    # read blocks of particle data
    bin_quantity = read_block_by_name(filename, blockname,
                                      info=info[getfield.(info, :block_name) .== blockname][1],
          					  		  parttype=0)

    x = read_block_by_name(filename, "POS",
			         	   info=info[getfield.(info, :block_name) .== "POS"][1],
				   		   parttype=0)

    rho = read_block_by_name(filename, "RHO",
    				   info=info[getfield.(info, :block_name) .== "RHO"][1],
    				   parttype=0)

    hsml = read_block_by_name(filename, "HSML",
    				      info=info[getfield.(info, :block_name) .== "HSML"][1],
    					parttype=0)

    m = read_block_by_name(filename, "MASS", parttype=0)

    if kernel_name == "WC6"
        kernel = WendlandC6()
    elseif kernel_name == "WC4"
	  kernel = WendlandC4()
    elseif kernel_name == "Quntic"
	  kernel = Quintic()
    elseif kernel_name == "Cubic"
	  kernel = Cubic()
    end

    if center_pos == [123456.7, 123456.7, 123456.7]
		if verbose
        	@info "Calculating COM..."
    	end
        center_pos = calculate_center_of_mass(x, m)
		println("COM: $center_pos")
    end

    if [dx, dy, dz] == [0.0, 0.0, 0.0]
        dx = h.boxsize
	  	dy = h.boxsize
	  	dz = h.boxsize
    end

    max_size = maximum([dx, dy, dz])

    par = mappingParameters(x_lim=[center_pos[1] - dx/2.0, center_pos[1] + dx/2.0],
    				    y_lim=[center_pos[2] - dx/2.0, center_pos[2] + dx/2.0],
    				    z_lim=[center_pos[3] - dx/2.0, center_pos[3] + dx/2.0],
    				    pixelSideLength=(max_size/2.0))

    if verbose
    	@info "Initial compilation run..."
    end
    d = sphMapping_2D(x, hsml, m, rho, bin_quantity,
					  param=par, kernel=kernel,
					  conserve_quantities=conserve_quantities,
					  show_progress=false)

    par = mappingParameters(x_lim=[center_pos[1] - dx/2.0, center_pos[1] + dx/2.0],
    				    y_lim=[center_pos[2] - dx/2.0, center_pos[2] + dx/2.0],
    				    z_lim=[center_pos[3] - dx/2.0, center_pos[3] + dx/2.0],
    				    pixelSideLength=(max_size/resolution))

    if verbose
		@info "Mapping..."
    end
    d = sphMapping_2D(x, hsml, m, rho, bin_quantity,
					  param=par, kernel=kernel,
					  conserve_quantities=conserve_quantities,
					  show_progress=true)

    if plot
		GR.imshow(d)
    end

    return d

end



"""

"""

function sphMapping(Pos, HSML, M, ρ, Bin_Quant;
                    param::mappingParameters, kernel::SPHKernel,
                    show_progress::Bool=true,
					conserve_quantities::Bool=true,
					dimensions::Int=2)

	if (dimensions == 2)
		return sphMapping_2D(Pos, HSML, M, ρ, Bin_Quant;
		                     param=param, kernel=kernel,
							 conserve_quantities=conserve_quantities,
		                     show_progress=show_progress)
	elseif (dimensions == 3 )
		return sphMapping_3D(Pos, HSML, M, ρ, Bin_Quant;
		                     param=param, kernel=kernel,
		                     show_progress=show_progress)
	end



end

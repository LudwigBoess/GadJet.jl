import GR
using Base.Threads
using Distributed

@everywhere include("mapping_functions.jl")

"""
    glimpse(filename::String, blockname::String[, ... ])

Reads relevant data from snapshot file and maps the quantity in block `blockname`
to a grid.

# Arguments
- `filename::String`: Name of snapshot file.
- `blockname::String`: Name of block that should be mapped.
- `center_pos::Array{Float64,1}=[123456.7, 123456.7, 123456.7]`: Image center. If not given center of mass is calculated.
- `dx::Float64`= 0.0`: Extent in x-direction. If `dx = 0.0` whole box is mapped.
- `dy::Float64`= 0.0`: Extent in y-direction. If `dy = 0.0` whole box is mapped.
- `dz::Float64`= 0.0`: Extent in z-direction. If `dz = 0.0` whole box is mapped.
- `kernel_name::String="WC6"`: Which kernel should be used. ["Cubic", "Quintic", "WC4", "WC6"]
- `resolution::Int64=500`: Number of pixels in the longest dimension.
- `run_dummy::Bool=true`: If a compilation run with 4 pixels should be performed.
- `parallel::Bool=true`: Run on multiple processors.
- `conserve_quantities::Bool=true`: If quantities should be conserved while mapping, like in Smac (Dolag et. al. 2005).
- `verbose::Bool=true`: Output information to console and show progress bar.
- `plot::Bool=false`: Plot the resulting map with GR.imshow()
"""
function glimpse(filename::String, blockname::String,
                 center_pos::Array{Float64,1}=[123456.7, 123456.7, 123456.7],
                 dx::Float64=0.0, dy::Float64=0.0, dz::Float64=0.0;
			     kernel_name::String="WC6",
                 resolution::Int64=500, run_dummy::Bool=true,
				 parallel::Bool=true,
				 conserve_quantities::Bool=false,
			     verbose::Bool=true, plot::Bool=false)


    if verbose
		@info "Reading data..."
    end

	filebase = filename

	# if the snapshot does not exist it may be split into multiple files
    if !isfile(filebase)

        if verbose
            @info "File: $filebase not found, looking for sub-files."
        end

        # try reading the first of the distributed snapshots
        filename = filebase * ".0"

        # throw error if file does not exist
        if !isfile(filename)
            error("File: $filename not found!")
        end

        h = head_to_obj(filename)

        if verbose
            @info "$(h.num_files) sub-files found."
        end

        nfiles = h.num_files
    else
        nfiles = 1
    end

    # read header of snapshot
    h = head_to_obj(filename)
    info = read_info(filename, verbose=false)

    if info == 1
        error("Glimpse only works for snapshots with info block, sorry!")
    end

	if nfiles == 1
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
	else

		if center_pos == [123456.7, 123456.7, 123456.7]
			error("glimpse for multiple files works only with a given central position!")
		end

		if dx == 0.0 || dy == 0.0 || dz == 0.0
			error("glimpse for multiple files works only with a given extent!")
		end

		blocks = [blockname, "POS", "RHO", "HSML", "MASS"]
		unique!(blocks)

		x0 = center_pos - 0.5 .* [dx, dy, dz]
		x1 = center_pos + 0.5 .* [dx, dy, dz]

		d = read_particles_in_box(filebase, blocks, x0, x1, verbose=verbose)

		bin_quantity = d[blockname]
		x 	 = d["POS"]
		rho  = d["RHO"]
		hsml = d["HSML"]
		m 	 = d["MASS"]

		# "deallocate" d
		d = nothing
	end

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

    if run_dummy
	    par = mappingParameters(center = center_pos,
								x_size = dx, y_size = dy, z_size = dz,
	    				    	Npixels = 2)

	    if verbose
	    	@info "Initial compilation run..."
	    end
	    d = sphMapping(x, hsml, m, rho, bin_quantity,
						  param=par, kernel=kernel,
						  conserve_quantities=conserve_quantities,
						  show_progress=false)
	end

    par = mappingParameters(center = center_pos,
							x_size = dx, y_size = dy, z_size = dz,
							Npixels = resolution)

    if verbose
		@info "Mapping..."
    end
    d = sphMapping(x, hsml, m, rho, bin_quantity,
					  param=par, kernel=kernel,
					  conserve_quantities=conserve_quantities,
					  show_progress=true)

    if plot
		GR.imshow(d)
    end

    return d

end


"""
    sphMapping(Pos, HSML, M, ρ, Bin_Quant;
	           param::mappingParameters,
			   kernel::SPHKernel[,
			   show_progress::Bool=true,
			   conserve_quantities::Bool=true,
			   parallel::Bool=true,
			   dimensions::Int=2])

Maps the data in `Bin_Quant` to a grid. Parameters of mapping are supplied in
`param` and the kernel to be used in `kernel`.

# Arguments
- `Pos`: Array with particle positions.
- `HSML`: Array with particle hsml.
- `M`: Array with particle masses.
- `ρ`: Array with particle densities.
- `Bin_Quant`: Array with particle quantity to be mapped.
- `kernel::SPHKernel`: Kernel object to be used.
- `show_progress::Bool=true`: Show progress bar.
- `parallel::Bool=true`: Run on multiple processors.
- `conserve_quantities::Bool=true`: If quantities should be conserved while mapping, like in Smac (Dolag et. al. 2005).
- `dimensions::Int=2`: Number of mapping dimensions (2 = to grid, 3 = to cube).
"""
function sphMapping(Pos, HSML, M, ρ, Bin_Quant;
                    param::mappingParameters,
					kernel::SPHKernel,
                    show_progress::Bool=true,
					conserve_quantities::Bool=true,
					parallel::Bool=true,
					dimensions::Int=2)

	if (dimensions == 2)

		# if !parrallel
		 	return sphMapping_2D(Pos, HSML, M, ρ, Bin_Quant;
		 	                     param=param, kernel=kernel,
								 conserve_quantities=conserve_quantities,
		 	                     show_progress=show_progress)
		# else
		#
		# 	N = length(M)
		# 	futures = Array{Future}(undef, nworkers())
		#
		# 	batches = Array{typeof(1:2)}(undef, nworkers())
		#
		# 	size = Int(floor(N/nworkers()))
		#
		# 	@inbounds for i = 1:nworkers()-1
		# 	    batches[i] = 1+(i-1)*size:i*size
		# 	end
		# 	batches[nworkers()] = 1+(nworkers()-1)*size:N
		#
		# 	# start remote processes
		# 	for (i, id) in enumerate(workers())
		# 		futures[i] = @spawnat id sphMapping_2D(Pos[batch[i],:], HSML[batch[i]],
		# 											   M[batch[i]], ρ[batch[i]],
		# 											   Bin_Quant[batch[i]];
		# 								   			   param=param, kernel=kernel,
		# 								   			   conserve_quantities=conserve_quantities,
		# 								   			   show_progress=false)
		# 	end
		#
		# 	return sum(fetch.(futures))
		# end

	elseif (dimensions == 3 )
		return sphMapping_3D(Pos, HSML, M, ρ, Bin_Quant;
		                     param=param, kernel=kernel,
		                     show_progress=show_progress)
	end



end

__precompile__()

module GadJet

    # general utility stuff
    include(joinpath(dirname(@__FILE__), "utility", "gravity_utility.jl"))

    # functions to read snapshots
    include(joinpath(dirname(@__FILE__), "read_snapshot", "gadget_types.jl"))
    #include(joinpath(dirname(@__FILE__), "read_snapshot", "dict_functions.jl"))
    #include(joinpath(dirname(@__FILE__), "read_snapshot", "obj_functions.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_header.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "snapshot_utilities.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_format_1.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_format_2.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_snapshot.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_subfind.jl"))
    include(joinpath(dirname(@__FILE__), "read_snapshot", "read_particle_in_box.jl"))

    # functions to write snapshots
    include(joinpath(dirname(@__FILE__), "write_snapshot", "write_snap.jl"))

    # unit conversion
    include(joinpath(dirname(@__FILE__), "unit_conversion", "unit_types.jl"))

    # sph to grid mapping internal module
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "kernels.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_types.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "mapping_functions.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_to_grid.jl"))

    # sph to grid mapping with Smac
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "smac1_utility.jl"))
    # sph to grid mapping with P-Smac2
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "smac2_utility.jl"))

    # riemann solvers
    include(joinpath(dirname(@__FILE__), "ideal_solutions", "cr_dsa_models.jl"))
    include(joinpath(dirname(@__FILE__), "ideal_solutions", "cr_sod_shock_main.jl"))
    include(joinpath(dirname(@__FILE__), "ideal_solutions", "setup_riemann_parameters.jl"))
    #include(joinpath(dirname(@__FILE__), "ideal_solutions", "cr_sod_shock.jl"))

    # sedov solution
    include(joinpath(dirname(@__FILE__), "ideal_solutions", "sedov_solution.jl"))

    include(joinpath(dirname(@__FILE__), "bp_cr_utility", "cr_datatypes.jl"))
    include(joinpath(dirname(@__FILE__), "bp_cr_utility", "analysis_functions.jl"))
    include(joinpath(dirname(@__FILE__), "bp_cr_utility", "get_detailled_data.jl"))


    export Header, Info_Line,       # types
           head_to_dict,
           snap_to_dict,
           head_to_obj,
           print_blocks,
           read_info,
           block_present,
           read_snap,
           read_block_by_name,      # similar to readnew.pro by Klaus Dolag
           read_header,
           read_particles_in_box,
           read_particles_in_volume,

           # subfind read
           read_subfind_header,
           read_subfind,
           find_most_massive_halo,

           # write snapshot functions
           write_header,
           write_block,

           # utility stuff
           calculate_center_of_mass,

           # Kernels
           Cubic,
           Quintic,
           WendlandC4,
           WendlandC6,
           kernel_value_2D,
           kernel_value_3D,

           # internal sph mapping
           mappingParameters,
           sphMapping,
           sphMapping_2D,
           sphMapping_3D,
           glimpse,

           # helper functions and datatypes for Smac
           Smac1ImageInfo,
           read_smac1_binary_image,
           read_smac1_binary_info,
           write_smac1_par,
           # helper function for P-Smac2
           write_smac2_par,

           # unit conversion
           GadgetPhysicalUnits,
           @u_str,
           strip_unit,
           # old riemann solver
           #RiemannParameters,   # datatype for riemann parameters
           #RiemannSolution,     # datatype for riemann solution
           #solveHydroShock,      # function that solves a standard sod shock

           # Test for different riemann solvers
           RiemannParameters,    # helper function to set up solution
           solve,                # overloaded function to solve riemann problems
           find_xs_first_guess,  # helper function to find initial guess for shock compression

           get_sedov_solution,   # wrapper function to get sedov data and ideal solution from snapshot

           # datatypes and helper functions for BP_REAL_CRs
           CRShockData,          # datatype to analyse single shocked particle
           readSingleCRShockDataFromOutputFile, # as the name says
           CRMomentumDistributionConfig, # config parameters for momentum distribution function
           CRMomentumDistribution,
           getCRMomentumDistributionFromPartID, # function to get distribution function
           calculateCREnergyInCGS,
           calculateCRNumber,
           get_detailled_shock_data,
           get_detailled_Dpp_data,
           get_detailled_radiative_data,
           get_detailled_adiabatic_data

end

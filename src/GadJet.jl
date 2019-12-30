#__precompile__()

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
    # functions to write snapshots
    include(joinpath(dirname(@__FILE__), "write_snapshot", "write_snap.jl"))
    # sph to grid mapping internal module
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "kernels.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_types.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "mapping_functions.jl"))
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "sph_to_grid.jl"))

    # sph to grid mapping with Smac
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "smac1_utility.jl"))
    # sph to grid mapping with P-Smac2
    include(joinpath(dirname(@__FILE__), "sph_to_grid", "smac2_utility.jl"))
    # unit conversion
    include(joinpath(dirname(@__FILE__), "unit_conversion", "unit_types.jl"))
    # old riemann solver
    #include(joinpath(dirname(@__FILE__), "riemann_solvers", "riemann_solver.jl"))

    # test for new riemann solvers
    include(joinpath(dirname(@__FILE__), "riemann_solvers", "cr_dsa_models.jl"))
    include(joinpath(dirname(@__FILE__), "riemann_solvers", "cr_sod_shock_main.jl"))
    include(joinpath(dirname(@__FILE__), "riemann_solvers", "setup_riemann_parameters.jl"))
    #include(joinpath(dirname(@__FILE__), "riemann_solvers", "cr_sod_shock.jl"))

    include(joinpath(dirname(@__FILE__), "bp_cr_utility", "cr_datatypes.jl"))
    include(joinpath(dirname(@__FILE__), "bp_cr_utility", "analysis_functions.jl"))


    export Header, Info_Line,       # types
           head_to_dict,
           snap_to_dict,
           head_to_obj,
           print_blocks,
           read_info,
           read_snap,
           read_block_by_name,      # similar to readnew.pro by Klaus Dolag
           read_header,

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
           # internal sph mapping
           mappingParameters,
           sphAdaptiveMapping,
           sphCenterMapping,
           glimpse,
           
           # helper functions and datatypes for Smac
           Smac1ImageInfo,
           read_smac1_binary_image,
           read_smac1_binary_info,
           # helper function for P-Smac2
           write_smac2_par,

           GadgetUnitFactors,
           # old riemann solver
           #RiemannParameters,   # datatype for riemann parameters
           #RiemannSolution,     # datatype for riemann solution
           #solveHydroShock,      # function that solves a standard sod shock

           # Test for different riemann solvers
           RiemannParameters,    # helper function to set up solution
           solve,                # overloaded function to solve riemann problems

           # datatypes and helper functions for BP_REAL_CRs
           CRShockData,          # datatype to analyse single shocked particle
           readSingleCRShockDataFromOutputFile, # as the name says
           CRMomentumDistributionConfig, # config parameters for momentum distribution function
           CRMomentumDistribution,
           getCRMomentumDistributionFromPartID, # function to get distribution function
           calculateCREnergyInCGS,
           calculateCRNumber



end

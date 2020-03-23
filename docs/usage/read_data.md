Read Snapshot Data
=========


Reading the header
------------------

Reading the header block of the simulation can be done by using:

```julia
 h = read_header(filename::String)
```

Where `h` is the returned header object:

```julia
mutable struct Header
    npart::Vector{Int32}                # an array of particle numbers per type in this snapshot
    massarr::Vector{Float64}            # an array of particle masses per type in this snapshot - if zero: MASS block present
    time::Float64                       # time / scale factor of the simulation
    z::Float64                          # redshift of the simulation
    flag_sfr::Int32                     # 1 if simulation was run with star formation, else 0
    flag_feedback::Int32                # 1 if simulation was run with stellar feedback, else 0
    nall::Vector{UInt32}                # total number of particles in the simulation
    flag_cooling::Int32                 # 1 if simulation was run with cooling, else 0
    num_files::Int32                    # number of snapshots over which the simulation is distributed
    boxsize::Float64                    # total size of the simulation box
    omega_0::Float64                    # Omega matter
    omega_l::Float64                    # Omega dark enery
    h0::Float64                         # little h
    flag_stellarage::Int32              # 1 if simulation was run with stellar age, else 0
    flag_metals::Int32                  # 1 if simulation was run with metals, else 0
    npartTotalHighWord::Vector{UInt32}  # weird
    flag_entropy_instead_u::Int32       # 1 if snapshot U field contains entropy instead of internal energy, else 0
    flag_doubleprecision::Int32         # 1 if snapshot is in double precision, else 0
    flag_ic_info::Int32                 
    lpt_scalingfactor::Float32
    fill::Vector{Int32}                 # the HEAD block needs to be filled with zeros to have a size of 256 bytes
end
```

This is equivalent to:

```julia
 h = head_to_obj(filename::String)
```

If you want to read the header information into a dictionary you can use:

```julia
 h = head_to_dict(filename::String)
```

Reading a snapshot
------------------

### Full snapshot
If you want to read a simulation snapshot into memory with GadJet.jl, it's as easy as this:

```julia
    data = read_snap(filename)
```

This will return a dictionary with the header information in `data["Header"]` and the blocks sorted by particle type.

As an example, this is how you would access the positions of the gas particles:

```julia
    data["Parttype0"]["POS"]
```


### Specific blocks

Reading specific blocks only works with Format 2 at the moment.

If you only want to read a specific block for a single particle type, e.g. positions of gas particles, you can use the function with a specified blockname and particle type like so:

```julia
    pos = read_snap(filename, "POS", parttype=0)
```

This will return an array of the datatype of your simulation, usually Float32.

If the snapshot has no info block this will fail unfortunately.

You can still read the specific block by supplying a hand-constructed `Info_Line` object:

```julia
mutable struct Info_Line
    block_name::String              # name of the data block, e.g. "POS"
    data_type::DataType             # datatype of the block, e.g. Float32 for single precision, Float64 for double
    n_dim::Int32                    # number of dimensions of the block, usually 1 or 3
    is_present::Vector{Int32}       # array of flags for which particle type this block is present,
                                    # e.g. gas only:  [ 1, 0, 0, 0, 0, 0 ]
                                    # e.g. gas + BHs: [ 1, 0, 0, 0, 0, 1 ]
end
```

and passing that to the function:

```julia
    pos = read_block_by_name(filename, "POS", info=pos_info, parttype=0)
```

where `pos_info` is a `Info_Line` object.

I will collect some example `Info_Line` objects in a later release to be able to read some common blocks even without a info block.

### Getting snapshot infos

If you have a Format 2 snapshot and just want to know what blocks the snapshot contains you can use the function

```julia
    print_blocks(filename)
```

to get an output of all block names.

If your simulation contains an INFO block you can read the info lines into `Info_Line` object like so:

```julia
    info = read_info(filename, verbose=true)
```

This will return an Array of `Info_Line` objects. The optional keyword `verbose` additionally gives the same functionality as `print_blocks` and prints the names to the console.


Large Simulations
-----------------

For large simulations Gadget distributes snapshots over multiple files. These files contain particles associated with specific peano-hilber keys.

To get all particles within a subvolume of the simulation you can use the functions `read_particles_in_box(...)` or `read_particles_in_volume(...)`.

`read_particles_in_box(...)` takes a box defined by a lower-left corner and an upper-right corner, constructs the peano hilbert keys, selects the relevant files and reads the particles from these files into a dictionary.

```julia
function read_particles_in_box(filename::String, blocks::Vector{String},
                               corner_lowerleft,
                               corner_upperright;
                               parttype::Int=0,
                               verbose::Bool=true)

            (...)

end
```

You can define an array of blocks you want to read, these will be read in parallel with simple multi-threading.

`read_particles_in_volume(...)` is a simple wrapper around `read_particles_in_box(...)`, where you can define a central position and a radius around it and it will construct the box containing that sphere for you and read all particles in it.

```julia
function read_particles_in_volume(filename::String, blocks::Vector{String},
                                  center_pos::Vector{AbstractFloat},
                                  radius::AbstractFloat;
                                  parttype::Int=0,
                                  verbose::Bool=true)

            (...)
end
```

In both functions `parttype` defines the particle type to be read, as in the previous read functions and `verbose` gives console output.

### Filename

With the snapshots being distributed over multiple filenames you need to be careful with that keyword. In this case filename refers to the base-name. Assuming you want to read snapshot 140, which is in the snapshot directory 140 the filename is

```julia
filename = "path/to/your/snapshot/directories/snapdir_140/snap_140"
```

GadJet will then automatically loop through the sub-snapshots which end in ".0", ".1", ... , ".N".

### Example

If you want to, e.g. read positions, velocities, masses, density and hsml for all gas particles within the virial radius of the most massive halo of a simulation you can do this as follows.

Assuming `pos_halo` is the position of the center of mass of the halo and `r_vir` is its virial radius you read the data with

```julia
blocks = ["POS", "VEL", "MASS", "RHO", "HSML"]

data   = read_particles_in_volume(filename, blocks, pos_halo, r_vir,
                                  parttype=0,
                                  verbose=true)
```

This will return a dictionary with the blocks as keys and containing the arrays for the particles.

```julia
data["POS"]  # array of positions
data["RHO"]  # array of densities
(...)
```


Read Subfind Data
=========

Reading the header
------------------

Similarly to the normal snapshot you can read the header of the subfind output into a `SubfindHeader` object

```julia
struct SubfindHeader
    nhalos::Int32                       # number of halos in the output file
    nsubhalos::Int32                    # number of subhalos in the output file
    nfof::Int32                         # number of particles in the FoF
    ngroups::Int32                      # number of large groups in the output file
    time::Float64                       # time / scale factor of the simulation
    z::Float64                          # redshift of the simulation
    tothalos::UInt32                    # total number of halos over all output files
    totsubhalos::UInt32                 # total number of subhalos over all output files
    totfof::UInt32                      # total number of particles in the FoF
    totgroups::UInt32                   # total number of large groups over all output files
    num_colors::Int32                   # number of colors
    boxsize::Float64                    # total size of the simulation box
    omega_0::Float64                    # Omega matter
    omega_l::Float64                    # Omega dark enery
    h0::Float64                         # little h
    flag_doubleprecision::Int32         # 1 if snapshot is in double precision, else 0
    flag_ic_info::Int32
end
```

using

```julia
h = read_subfind_header(filename::String)
```

Reading the subfind files
-------------------------

If you compiled Gadget with `WRITE_SUB_IN_SNAP_FORMAT` you can read the subfind output like you would a normal snapshot, with any of the above methods. For convenience you can also use a helper function provided by GadJet. Since each of the blocks is only relevant for either halos, subhalos, Fof or large groups you don't need to define a particly type, aka halo type in this case.

So in order to read the virial radius of the halos in a file you can simply use

```julia
R_vir = read_subfind(filename, "RVIR")
```

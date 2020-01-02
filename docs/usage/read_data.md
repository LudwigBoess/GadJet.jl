Read Data
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

"""
    Functions in this file read the subfind output.

"""

struct SubfindHeader
    nhalos::Int32                # an array of particle numbers per type in this snapshot
    nsubhalos::Int32
    nfof::Int32
    ngroups::Int32
    time::Float64                       # time / scale factor of the simulation
    z::Float64                          # redshift of the simulation
    tothalos::UInt32
    totsubhalos::UInt32
    totfof::UInt32
    totgroups::UInt32
    num_colors::Int32                   # number of colors
    boxsize::Float64                    # total size of the simulation box
    omega_0::Float64                    # Omega matter
    omega_l::Float64                    # Omega dark enery
    h0::Float64                         # little h
    flag_doubleprecision::Int32         # 1 if snapshot is in double precision, else 0
    flag_ic_info::Int32
end

function read_subfind_header(filename::String)

    f = open(filename)
    blocksize = read(f, Int32)

    if blocksize == 8
        swap = 0
        snap_format = 2
    elseif blocksize == 256
        swap = 0
        snap_format = 1
    else
        blocksize = bswap(blocksize)
        if blocksize == 8
            swap = 1
            snap_format = 2
        elseif blocksize == 256
            swap = 1
            snap_format = 1
        else
            println("incorrect file format encountered when reading header of", filename)
        end
    end

    #println("Reading snapshot format: ", snap_format)

    if snap_format == 2
        seek(f, 16)
        skip_line = read(f, Int32)
    end

    nhalos = read(f, Int32)
    nsubhalos = read(f, Int32)
    nfof = read(f, Int32)
    ngroups = read(f, Int32)

    # rest of numbers and massarr are obsolete
    DUMMY = read!(f, Array{Int32,1}(undef,2))
    DUMMY = read!(f, Array{Float64,1}(undef,6))

    time = read(f, Float64)
    z = read(f, Float64)

    # flags not needed
    DUMMY = read(f, Int32)
    DUMMY = read(f, Int32)

    tothalos = read(f, UInt32)
    totsubhalos = read(f, UInt32)
    totfof = read(f, UInt32)
    totgroups = read(f, UInt32)

    DUMMY = read!(f, Array{Int32,1}(undef,2))
    DUMMY = read(f, Int32)

    num_files = read(f, Int32)
    boxsize = read(f, Float64)
    omega_0 = read(f, Float64)
    omega_l = read(f, Float64)
    h0 = read(f, Float64)

    DUMMY = read(f, Int32)
    DUMMY = read(f, Int32)
    DUMMY = read!(f, Array{UInt32,1}(undef,6))
    DUMMY = read(f, Int32)

    flag_doubleprecision = read(f, Int32)
    flag_ic_info = read(f, Int32)

    close(f)

    return SubfindHeader(nhalos, nsubhalos, nfof, ngroups,
                         time, z,
                         tothalos, totsubhalos, totfof, totgroups,
                         num_files,
                         boxsize, omega_0, omega_l, h0,
                         flag_doubleprecision,
                         flag_ic_info)
end


function read_subfind(filename::String, blockname::String)

    info = read_info(filename)

    info = info[getfield.(info, :block_name) .== blockname][1]

    parttype = findall(info.is_present .== 1)[1] - 1

    return read_snap(filename, blockname, parttype)

end

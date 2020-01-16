"""
        Helper function for Smac reading of binary images
"""

struct Smac1ImageInfo

    snap::Int32                 # number of input snapshot
    z::Float32                  # redshift of snapshot
    m_vir::Float32              # virial mass of halo
    r_vir::Float32              # virial radius of halo
    xcm::Float32                # x coordinate of image center
    ycm::Float32                # y coordinate of image center
    zcm::Float32                # z coordinate of image center
    z_slice_kpc::Float32        # depth of the image in kpc
    boxsize_kpc::Float32        # xy-size of the image in kpc
    boxsize_pix::Float32        # xy-size of the image in pixels
    pixsize_kpc::Float32        # size of one pixel in kpc
    xlim::Array{Float64,1}      # x limits of image
    ylim::Array{Float64,1}      # y limits of image
    zlim::Array{Float64,1}      # z limits of image
    units::String               # unitstring of image

    function Smac1ImageInfo(snap::Int32, z::Float32, m_vir::Float32, r_vir::Float32,
                       xcm::Float32, ycm::Float32, zcm::Float32,
                       z_slice_kpc::Float32,
                       boxsize_kpc::Float32, boxsize_pix::Float32, pixsize_kpc::Float32,
                       units::String)

        xlim = [xcm - boxsize_kpc/2.0, xcm + boxsize_kpc/2.0]
        ylim = [ycm - boxsize_kpc/2.0, ycm + boxsize_kpc/2.0]
        zlim = [zcm - z_slice_kpc/2.0, zcm + z_slice_kpc/2.0]

        new(snap, z, m_vir, r_vir,
            xcm, ycm, zcm,
            z_slice_kpc,
            boxsize_kpc, boxsize_pix, pixsize_kpc,
            xlim, ylim, zlim, units)
    end

end

function read_smac1_binary_image(fi)
    f = open(fi)
    first = read(f, Int32)
    n_pixels = Int64(sqrt.(first/4.0))
    a = read!(f, Array{Float32,2}(undef,n_pixels,n_pixels))
    last = read(f, Int32)
    close(f)
    if first == last
        return collect(transpose(a))
    else
        error("Read error: Incorrect image format!")
    end
end

function read_smac1_binary_info(fi)

    f = open(fi)

    # skip image
    image = read(f, Int32)
        seek(f, position(f) + image + 8)

    # read info block
    snap = read(f, Int32)
        seek(f, position(f) + 8)
    z = read(f, Float32)
        seek(f, position(f) + 8)
    m_vir = read(f, Float32)
        seek(f, position(f) + 8)
    r_vir = read(f, Float32)
        seek(f, position(f) + 8)
    xcm = read(f, Float32)
        seek(f, position(f) + 8)
    ycm = read(f, Float32)
        seek(f, position(f) + 8)
    zcm = read(f, Float32)
        seek(f, position(f) + 8)
    z_slice_kpc = read(f, Float32)
        seek(f, position(f) + 8)
    boxsize_kpc = read(f, Float32)
        seek(f, position(f) + 8)
    boxsize_pix = read(f, Float32)
        seek(f, position(f) + 8)
    pixsize_kpc = read(f, Float32)
        seek(f, position(f) + 4)
    # read length of unit string
    unitstring_length = read(f, Int32)
    # read unti string
    unit = String(Char.((read!(f, Array{Int8,1}(undef,unitstring_length)))))
    close(f)

    img = Smac1ImageInfo(snap, z, m_vir, r_vir,
                    xcm, ycm, zcm,
                    z_slice_kpc,
                    boxsize_kpc, boxsize_pix, pixsize_kpc,
                    unit)

    return img
end

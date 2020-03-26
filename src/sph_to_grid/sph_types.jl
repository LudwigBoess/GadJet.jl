"""
            Objects for SPH mapping.


    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""


mutable struct mappingParameters

    x_lim::Vector{Float64}
    y_lim::Vector{Float64}
    z_lim::Vector{Float64}
    x_max::Vector{Float64}
    y_max::Vector{Float64}
    z_max::Vector{Float64}
    pixelSideLength::Float64
    pixelArea::Float64
    Npixels::Int64
    xy_size::Float64
    z_size::Float64
    x::Array{Float64,1}
    y::Array{Float64,1}
    z::Array{Float64,1}

    function mappingParameters(;x_lim::Vector{Float64},
                                 y_lim::Vector{Float64},
                                 z_lim::Vector{Float64},
                                 pixelSideLength::Float64=0.0,
                                 Npixels::Int64=0)

        xy_size = maximum( [ abs(x_lim[1] - x_lim[2]),
                             abs(y_lim[1] - y_lim[2]) ] )

        z_size = abs(z_lim[1]) + abs(z_lim[2])

        if (pixelSideLength == 0.0) & (Npixels != 0)
            pixelSideLength = xy_size/Npixels
        elseif (pixelSideLength != 0.0) & (Npixels == 0)
            Npixels = floor(Int64, xy_size/pixelSideLength)
        else
            error("Please specify pixelSideLenght or number of pixels!")
        end

        pixelArea = pixelSideLength^2

        x_max = collect(x_lim[1]:pixelSideLength:x_lim[2])
        y_max = collect(y_lim[1]:pixelSideLength:y_lim[2])
        z_max = collect(z_lim[1]:pixelSideLength:z_lim[2])

        x = collect(x_lim[1] : pixelSideLength/2. : x_lim[2])[2:2:end]
        y = collect(y_lim[1] : pixelSideLength/2. : y_lim[2])[2:2:end]
        z = collect(z_lim[1] : pixelSideLength/2. : z_lim[2])[2:2:end]


        new(x_lim, y_lim, z_lim,
            x_max, y_max, z_max,
            pixelSideLength,
            pixelArea,
            Npixels,
            xy_size, z_size,
            x, y, z)

    end

end

# par = mappingParameters2(x_lim=[-10.0 , 10.], y_lim=[-10.0 , 10.], z_lim=[-2.0 , 2.0],
#                          pixelSideLength=0.01, normalize=false, broadening=false)

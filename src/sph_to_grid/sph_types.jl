"""
            Objects for SPH mapping.


    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""

#using Distributed



"""
    mappingParameters(; x_lim::Vector{Float64},
                        y_lim::Vector{Float64},
                        z_lim::Vector{Float64},
                        pixelSideLength::Float64=0.0,
                        Npixels::Int64=0)

Parameter object for sph to grid mapping. Define either `pixelSideLength` of `Npixels`.
"""
struct mappingParameters

    x_lim::Vector{Float64}
    y_lim::Vector{Float64}
    z_lim::Vector{Float64}
    center::Vector{Float64}
    pixelSideLength::Float64
    pixelArea::Float64
    Npixels::Int64
    x_size::Float64
    y_size::Float64
    z_size::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}

    function mappingParameters(;x_lim::Vector{Float64}   = [-1.0, -1.0],
                                y_lim::Vector{Float64}   = [-1.0, -1.0],
                                z_lim::Vector{Float64}   = [-1.0, -1.0],
                                center::Vector{Float64}  = [-1.0, -1.0, -1.0],
                                x_size::Float64          =  -1.0,
                                y_size::Float64          =  -1.0,
                                z_size::Float64          =  -1.0,
                                pixelSideLength::Float64 =  -1.0,
                                Npixels::Int64           =   0)


        # calculate limits if center position and sizes are given
        if ( x_lim == [-1.0, -1.0] && ( y_lim == [-1.0, -1.0] && z_lim == [-1.0, -1.0] ))

            if (center != [-1.0, -1.0, -1.0]) && ( x_size != -1.0 && (y_size != -1.0 && z_size != -1.0) )

               x_lim = [ center[1] - 0.5x_size, center[1] + 0.5x_size ]
               y_lim = [ center[2] - 0.5y_size, center[2] + 0.5y_size ]
               z_lim = [ center[3] - 0.5z_size, center[3] + 0.5z_size ]

            else
                error("Giving a center position requires extent in x, y and z direction.")
            end
        end

        # calculate side lengths from limits
        if x_size == -1.0
            x_size = x_lim[2] - x_lim[1]
        end
        if y_size == -1.0
            y_size = y_lim[2] - y_lim[1]
        end
        if z_size == -1.0
            z_size = z_lim[2] - z_lim[1]
        end

        # find the maximum extent of the map
        max_size = max(x_size, y_size)

        if (pixelSideLength == -1.0) & (Npixels != 0)
            pixelSideLength = max_size/Npixels
        elseif (pixelSideLength != 0.0) & (Npixels == 0)
            Npixels = floor(Int64, max_size/pixelSideLength)
        else
            error("Please specify pixelSideLenght or number of pixels!")
        end

        pixelArea = pixelSideLength^2

        x = collect(range(x_lim[1]+0.5pixelSideLength, x_lim[2]-0.5pixelSideLength, step = pixelSideLength))
        y = collect(range(y_lim[1]+0.5pixelSideLength, y_lim[2]-0.5pixelSideLength, step = pixelSideLength))
        z = collect(range(z_lim[1]+0.5pixelSideLength, z_lim[2]-0.5pixelSideLength, step = pixelSideLength))


        new(x_lim, y_lim, z_lim,
            center,
            pixelSideLength,
            pixelArea,
            Npixels,
            x_size, y_size, z_size,
            x, y, z)

    end
end

# mutable struct mappingParameters
#
#     pixelsPerDimension::Int64
#     centerX::Float64
#     centerY::Float64
#     centerZ::Float64
#     minZ::Float64
#     maxZ::Float64
#     pixelSideLength::Float64
#     lineOfSight::String
#     method::String
#     preview::Bool
#     normalize::Bool
#     broadening::Bool
#     x::Array{Float64,1}
#     y::Array{Float64,1}
#     c::Array{Float64,2}
#
#
#     function mappingParameters(pixelsPerDimension::Int64,
#                                centerX::Float64, centerY::Float64, centerZ::Float64,
#                                minZ::Float64, maxZ::Float64,
#                                pixelSideLength::Float64,
#                                lineOfSight::String,
#                                method::String,
#                                preview::Bool=false,
#                                normalize::Bool=false,
#                                broadening::Bool=false)
#
#         new(pixelsPerDimension,
#             centerX, centerY, centerZ,
#             minZ, maxZ,
#             pixelSideLength,
#             lineOfSight,
#             method,
#             preview,
#             normalize,
#             broadening)
#
#     end
#
# end

mutable struct mappingParameters

    x_lim::Vector{Float64}
    y_lim::Vector{Float64}
    z_lim::Vector{Float64}
    x_max::Vector{Float64}
    y_max::Vector{Float64}
    z_max::Vector{Float64}
    pixelSideLength::Float64
    normalize::Bool
    broadening::Bool
    x::Array{Float64,1}
    y::Array{Float64,1}
    z::Array{Float64,1}

    function mappingParameters(;x_lim::Vector{Float64},
                                 y_lim::Vector{Float64},
                                 z_lim::Vector{Float64},
                                 pixelSideLength::Float64,
                                 normalize::Bool=false,
                                 broadening::Bool=true)

        x_max = collect(x_lim[1]:pixelSideLength:x_lim[2])
        y_max = collect(y_lim[1]:pixelSideLength:y_lim[2])
        z_max = collect(z_lim[1]:pixelSideLength:z_lim[2])

        x = collect(x_lim[1] : pixelSideLength/2. : x_lim[2])[2:2:end]
        y = collect(y_lim[1] : pixelSideLength/2. : y_lim[2])[2:2:end]
        z = collect(z_lim[1] : pixelSideLength/2. : z_lim[2])[2:2:end]


        new(x_lim, y_lim, z_lim,
            x_max, y_max, z_max, 
            pixelSideLength,
            normalize,
            broadening,
            x, y, z)



    end

end

# par = mappingParameters2(x_lim=[-10.0 , 10.], y_lim=[-10.0 , 10.], z_lim=[-2.0 , 2.0],
#                          pixelSideLength=0.01, normalize=false, broadening=false)

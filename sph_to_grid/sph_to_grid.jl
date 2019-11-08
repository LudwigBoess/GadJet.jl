"""
            Functions for sph mapping to grid.

            sphCenterMapping works like SPLASH.

            sphAdaptiveMapping is adapted from SPHMapper by Dr. Alexander Arth
            (private conversations). A similar implementation is used in Pygad.
            Publication: https://arxiv.org/abs/1803.03652

            !!!! Not done yet !!!!


    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""

include(joinpath(dirname(@__FILE__), "kernels.jl"))
include(joinpath(dirname(@__FILE__), "sph_types.jl"))

using Statistics
using ProgressMeter
using Base.Threads
#using SharedArrays

@inline function find_min_pixel(pos::Float64, hsml::Float64,
                                minCoord::Float64, pixSize::Float64)

    pixmin = floor((pos - hsml - minCoord) / pixSize )
    if pixmin < 1
        return 1
    else
        return pixmin
    end
end

@inline function find_max_pixel(pos::Float64, hsml::Float64,
                                minCoord::Float64, pixSize::Float64,
                                max::Int64)

    pixmax = floor((pos + hsml - minCoord) / pixSize )
    if pixmax > max
        return max
    else
        return pixmax
    end
end

@inline function check_in_image(pos::Float64, hsml::Float64,
                                minCoord::Float64, maxCoord::Float64)

    if ( (minCoord - hsml) <= pos <= (maxCoord + hsml) )
        return true
    else
        return false
    end

end

function sphCenterMapping(Pos::Array{Float64,2}, HSML::Array{Float64,2}, M::Array{Float64,2},
                          ρ::Array{Float64,2}, Bin_Quant::Array{Float64,2};
                          param::mappingParameters, kernel)

    N = length(M)  # number of particles

    #val = SharedArray{Float64,2}(length(param.x), length(param.y))
    val = zeros(length(param.x), length(param.y))

    minCoords = [param.x[1], param.y[1], param.z[1]]
    maxCoords = [param.x[end], param.y[end], param.z[end]]

    max_pixel = [length(param.x), length(param.y), length(param.z)]

    @threads for p = 1:N

        in_image = false

        @inbounds for dim = 1:3
            in_image = check_in_image(Pos[p,dim], HSML[p],
                                      minCoords[dim], maxCoords[dim])
        end

        if in_image

            pixmin = Vector{Int64}(undef,3)
            pixmax = Vector{Int64}(undef,3)

            @inbounds for dim = 1:3

                pixmin[dim] = find_min_pixel(Pos[p,dim], HSML[p], minCoords[dim],
                                             param.pixelSideLength)

                pixmax[dim] = find_max_pixel(Pos[p,dim], HSML[p], minCoords[dim],
                                             param.pixelSideLength, max_pixel[dim])
            end

            @inbounds for i = pixmin[1]:pixmax[1]
                @inbounds for j = pixmin[2]:pixmax[2]
                    @inbounds for k = pixmin[3]:pixmax[3]

                        distance = sqrt( (param.x[i] - Pos[p,1])^2 +
                                         (param.y[j] - Pos[p,2])^2 +
                                         (param.z[k] - Pos[p,3])^2 )

                        val[i,j] += Bin_Quant[p] * M[p] / ρ[p] * kernel_value(kernel, distance/HSML[p], HSML[p])

                    end # end z-loop
                end # end y-loop
            end # end x-loop

        end # end check if in image

    end

   # c = zeros(length(param.y), length(param.x))
   #
   # cell_count = 1
   # for i = 1:length(param.x)
   #     for j = 1:length(param.y)
   #
   #         c[i,j] = sum(val[i,j,:])
   #
   #     end
   # end

   return val
end

par = mappingParameters(x_lim=[0.0,6.0], y_lim=[0.0,6.0], z_lim=[2.9,3.1],
                        pixelSideLength=0.02)
kernel = WendlandC6(295)


function get_bounds(pos::Array{Float32,2}, min::Float64, max::Float64, axis::Int64,
                    range_arr::Vector{Int64})

    k = findall(pos[range_arr,axis] .>= min  )
    m = findall(pos[range_arr[k],axis] .<= max  )

    return range_arr[k[m]]

end


function sphAdaptiveMapping(Pos, HSML, M, ρ, Bin_Quant; param::mappingParameters, kernel)

        """ Integral conserving SPH mapping. https://arxiv.org/abs/1803.03652
            !!! Diverges at the moment for large particle numbers due to lack
            of kernel normalisation !!!
        """


    N = length(M)
    pix_Vol = param.pixelSideLength^3
    boundary = false

    # println("mapping ", N, " particles")
    # println("Running on ", Threads.nthreads(), " Threads.")


    val = zeros(length(param.x), length(param.y),length(param.z))



    minCoords = [param.x[1], param.y[1], param.z[1]]

    #@showprogress 1 "Mapping..." for p = 1:N
    @threads for p = 1:N

        @inbounds pos = Pos[p,:]
        @inbounds hsml = HSML[p]
        @inbounds mass = M[p]
        @inbounds rho = ρ[p]
        @inbounds bin_quant = Bin_Quant[p]

        pixmin = Vector{Int64}(undef,3)
        pixmax = Vector{Int64}(undef,3)

        for dim = 1:3

            @inbounds pixmin[dim] = Int64( round((pos[dim] - hsml - minCoords[dim]) / param.pixelSideLength ))

            if pixmin[dim] < 1

                boundary = true
                pixmin[dim] = 1

            end

            pixmax[dim] = Int64( round((pos[dim] + hsml - minCoords[dim]) / param.pixelSideLength ))

            if dim == 1
                max = length(param.x)
            elseif dim == 2
                max = length(param.y)
            else
                max = length(param.z)
            end


            if pixmax[dim] > max
                boundary = true
                pixmax[dim] = max
            end

        end

        norm = 0.
        loop = true

        if boundary && hsml > 0.5 * param.pixelSideLength
            norm = 1.0
        else
            # closestPixel = Vector{Int64}(undef,3)
            # closestDistance = Inf

            for i = pixmin[1]:pixmax[1]
                for j = pixmin[2]:pixmax[2]
                    for k = pixmin[3]:pixmax[3]

                        @inbounds distance = sqrt.( (param.x[i] - pos[1])^2 +
                                                    (param.y[j] - pos[2])^2 +
                                                    (param.z[k] - pos[3])^2 )

                        # if distance < closestDistance
                        #     closestPixel = [i, j, k]            # ?
                        #     closestDistance = distance          # ?
                        # end

                        norm += kernel_value(kernel, distance/hsml, Float64.(hsml))

                    end
                end
            end

            norm *= pix_Vol

            if norm < 0.001
                # Line 468 ???
                #val = mass / ( rho * pix_Vol)
            else
                norm = 1.0/norm
            end

        end

        if loop == true        # ??

            closestPixel = Vector{Int64}(undef,3)
            closestDistance = Inf

            for i = pixmin[1]:pixmax[1]
                for j = pixmin[2]:pixmax[2]
                    for k = pixmin[3]:pixmax[3]

                        @inbounds distance = sqrt.( (param.x[i] - pos[1])^2 +
                                          (param.y[j] - pos[2])^2 +
                                          (param.z[k] - pos[3])^2 )

                        if distance < closestDistance
                            closestPixel = [i, j, k]
                            closestDistance = distance
                        end

                        val[i,j,k] += bin_quant * mass / rho * kernel_value(kernel, distance/hsml, Float64.(hsml))  * norm

                    end
                end
            end

        end



    end
    #
    # x = zeros(length(param.x)*length(param.y))
    # y = zeros(length(param.x)*length(param.y))
    # c = zeros(length(param.x)*length(param.y))
    # x = zeros(length(param.x_max), length(param.y_max))
    # y = zeros(length(param.x_max), length(param.y_max))
    # c = zeros(length(param.x), length(param.y))
     c = zeros(length(param.y), length(param.x))

    #cell_count = 1
    for i = 1:length(param.x)
        for j = 1:length(param.y)


            @inbounds c[i,j] = sum(val[i,j,:])
            #c[j,i] = sum(val[i,j,:]) # transpose for imshow()

            # x[cell_count] = param.x[i]
            # y[cell_count] = param.y[j]
            # c[cell_count] = sum(val[i,j,:])
            #
            # cell_count += 1

        end
    end


    return c

end

# """
#     Mapping similar to that of of SPLASH by
# """
# function splash2Dmapping(Pos, HSML, M, ρ, Bin_Quant; param::mappingParameters, kernel)
# end

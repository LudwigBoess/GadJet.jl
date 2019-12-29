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

# include(joinpath(dirname(@__FILE__), "kernels.jl"))
# include(joinpath(dirname(@__FILE__), "sph_types.jl"))

using Statistics
using ProgressMeter
using Base.Threads
using SharedArrays

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

@inline function get_distance(Pos::Array{Float64,1}, x::Float64, y::Float64, z::Float64)

    return sqrt( (x - pos[1] )^2 + ( y - pos[2] )^2 + ( z - pos[3])^2 )
end


# function sphCenterMapping(Pos::Array{Float64,2}, HSML::Array{Float64,2}, M::Array{Float64,2},
#                           ρ::Array{Float64,2}, Bin_Quant::Array{Float64,2};
#                           param::mappingParameters, kernel)
#
#     N = length(M)  # number of particles
#
#     #val = SharedArray{Float64,2}(length(param.x), length(param.y))
#     println("test")
#     val = zeros(length(param.x), length(param.y))
#
#     minCoords = [param.x[1], param.y[1], param.z[1]]
#     maxCoords = [param.x[end], param.y[end], param.z[end]]
#
#     max_pixel = [length(param.x), length(param.y), length(param.z)]
#
#     @inbounds for p = 1:N
#
#         in_image = false
#
#         @inbounds for dim = 1:3
#
#             if ( (minCoords[dim] - HSML[p]) <= Pos[p,dim] <= (maxCoords[dim] + HSML[p]) )
#                 in_image = true
#             else
#                 in_image = false
#             end
#             # in_image = check_in_image(Pos[p,dim], HSML[p],
#             #                           minCoords[dim], maxCoords[dim])
#         end
#
#         if in_image
#
#             pixmin = Vector{Int64}(undef,3)
#             pixmax = Vector{Int64}(undef,3)
#
#             @inbounds for dim = 1:3
#
#                 pix = floor((Pos[p,dim] - HSML[p] - minCoords[dim]) / param.pixelSideLength )
#                 if pix < 1
#                     pixmin[dim] = 1
#                 else
#                     pixmin[dim] = pix
#                 end
#
#                 # pixmin[dim] = find_min_pixel(Pos[p,dim], HSML[p], minCoords[dim],
#                 #                              param.pixelSideLength)
#                 #
#
#                 pix = floor((Pos[p,dim] + HSML[p] - minCoords[dim]) / param.pixelSideLength )
#                 if pix > max_pixel[dim]
#                     pixmax[dim] = max_pixel[dim]
#                 else
#                     pixmax[dim] = pix
#                 end
#                 # pixmax[dim] = find_max_pixel(Pos[p,dim], HSML[p], minCoords[dim],
#                 #                              param.pixelSideLength, max_pixel[dim])
#             end
#
#             @inbounds for i = pixmin[1]:pixmax[1]
#                 @inbounds for j = pixmin[2]:pixmax[2]
#                     @inbounds for k = pixmin[3]:pixmax[3]
#
#                         distance = sqrt( (param.x[i] - Pos[p,1])^2 +
#                                          (param.y[j] - Pos[p,2])^2 +
#                                          (param.z[k] - Pos[p,3])^2 )
#
#                         val[i,j] += Bin_Quant[p] * M[p] / ρ[p] * kernel_value(kernel, distance/HSML[p], HSML[p])
#
#                     end # end z-loop
#                 end # end y-loop
#             end # end x-loop
#
#         end # end check if in image
#
#     end
#
#    return val
# end

@inline function get_d_hsml(dx::Float64, dy::Float64, dz::Float64, hsml::Float64)
    result::Float64 = sqrt( dx*dx + dy*dy + dz*dz ) / hsml
    return result
end

function sphCenterMapping(Pos::Array{Float64,2}, HSML::Array{Float64,2}, M::Array{Float64,2},
                          ρ::Array{Float64,2}, Bin_Quant::Array{Float64,2};
                          param::mappingParameters, kernel,
                          show_progress::Bool=true)

    N = length(M)  # number of particles

    val = SharedArray{Float64,2}(length(param.x), length(param.y))

    #val = zeros(Float64, length(param.x), length(param.y))

    minCoords = [param.x[1], param.y[1], param.z[1]]
    maxCoords = [param.x[end], param.y[end], param.z[end]]

    max_pixel = [length(param.x), length(param.y), length(param.z)]

    if show_progress
        P = Progress(N)
        idx = 0
        P_lock = SpinLock()
    end

    @threads for p = 1:N

        # save stuff from array to single variables
        pos     = Pos[p,:]
        hsml    = HSML[p]
        bin_q   = Bin_Quant[p]
        m       = M[p]
        rho     = ρ[p]


        in_image::Bool = false


        @inbounds for dim = 1:3

            if ( (minCoords[dim] - hsml) <= pos[dim] <= (maxCoords[dim] + hsml) )
                in_image = true
            else
                in_image = false
            end
            # in_image = check_in_image(Pos[p,dim], HSML[p],
            #                           minCoords[dim], maxCoords[dim])
        end

        if in_image

            pixmin = Vector{Int64}(undef,3)
            pixmax = Vector{Int64}(undef,3)

            @inbounds for dim = 1:3

                pix = floor(Int64, (pos[dim] - hsml - minCoords[dim]) / param.pixelSideLength )
                if pix < 1
                    pixmin[dim] = 1
                else
                    pixmin[dim] = pix
                end

                # pixmin[dim] = find_min_pixel(Pos[p,dim], HSML[p], minCoords[dim],
                #                              param.pixelSideLength)
                #

                pix = floor(Int64, (pos[dim] + hsml - minCoords[dim]) / param.pixelSideLength )
                if pix > max_pixel[dim]
                    pixmax[dim] = max_pixel[dim]
                else
                    pixmax[dim] = pix
                end
                # pixmax[dim] = find_max_pixel(Pos[p,dim], HSML[p], minCoords[dim],
                #                              param.pixelSideLength, max_pixel[dim])
            end

            bin_prefac::Float64 = bin_q * m / rho

            @inbounds for i = pixmin[1]:pixmax[1]
                dx::Float64 = param.x[i] - pos[1]
                @inbounds for j = pixmin[2]:pixmax[2]
                    dy::Float64 = param.y[j] - pos[2]
                    @inbounds for k = pixmin[3]:pixmax[3]
                        dz::Float64 = param.z[k] - pos[3]

                        distance_hsml::Float64 = get_d_hsml(dx, dy, dz, hsml)#sqrt( dx*dx + dy*dy + dz*dz ) / hsml

                        if distance_hsml < 1.0
                            val[i,j] += bin_prefac * kernel_value(kernel, distance_hsml, hsml)
                        end

                    end # end z-loop
                end # end y-loop
            end # end x-loop

        end # end check if in image

        if show_progress
            lock(P_lock)
            idx += 1
            ProgressMeter.update!(P, idx)
            unlock(P_lock)
        end
    end

   return val
end


function sphCenterMapping_toCube(Pos::Array{Float64,2}, HSML::Array{Float64,2},
                                 M::Array{Float64,2},
                                 ρ::Array{Float64,2}, Bin_Quant::Array{Float64,2};
                                 param::mappingParameters, kernel)

    N = length(M)  # number of particles

    val = zeros(length(param.x), length(param.y), length(param.z))

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

                        val[i,j,k] += Bin_Quant[p] * M[p] / ρ[p] * kernel_value(kernel, distance/HSML[p], HSML[p])

                    end # end z-loop
                end # end y-loop
            end # end x-loop

        end # end check if in image

    end

   return val
end

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

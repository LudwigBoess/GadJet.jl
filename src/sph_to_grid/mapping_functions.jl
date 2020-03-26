"""
            Functions for sph mapping to grid.

            sphCenterMapping works like SPLASH.

            sphAdaptiveMapping is adapted from SPHMapper by Dr. Alexander Arth
            (private conversations). A similar implementation is used in Pygad.
            Publication: https://arxiv.org/abs/1803.03652

    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""

using Statistics
using ProgressMeter
using Base.Threads
#using SharedArrays

@inline function get_d_hsml_2D(dx::Float64, dy::Float64, hsml_inv::Float64)
    sqrt( dx*dx + dy*dy ) * hsml_inv
end

@inline function get_d_hsml_3D(dx::Float64, dy::Float64, dz::Float64,
                               hsml_inv::Float64)
    sqrt( dx*dx + dy*dy + dz*dz ) * hsml_inv
end

@inline function check_in_image(pos::Float64, hsml::Float64, minCoords, maxCoords)

    if ( (minCoords - hsml) <= pos <= (maxCoords + hsml) )
        return true
    else
        return false
    end
end

@inline function find_min_pixel(pos::Float64, hsml::Float64,
                                minCoords::AbstractFloat,
                                pixelSideLength::Float64)

    pix = floor(Int64, (pos - hsml - minCoords) / pixelSideLength )

    return max(pix, 1)
end

@inline function find_max_pixel(pos::Float64, hsml::Float64,
                                minCoords::AbstractFloat, pixelSideLength::Float64,
                                max_pixel::Integer)

    pix = floor(Int64, (pos + hsml - minCoords) / pixelSideLength )

    return min(pix, max_pixel)
end


function sphMapping_2D(Pos, HSML, M, ρ, Bin_Quant;
                       param::mappingParameters, kernel::SPHKernel,
                       conserve_quantities::Bool=false,
                       show_progress::Bool=false)

    N = length(M)  # number of particles

    image = zeros(length(param.x), length(param.y))

    minCoords = [param.x[1], param.y[1], param.z[1]]
    maxCoords = [param.x[end], param.y[end], param.z[end]]

    max_pixel = [length(param.x), length(param.y)]

    particles_in_image = 0

    if show_progress
        P = Progress(N)
        idx = 0
        #P_lock = SpinLock()  # uncomment to make thread-safe if needed in the future
    end

    @inbounds for p = 1:N

        # save stuff from array to single variables
        @inbounds pos  = Float64.(Pos[p,:])
        @inbounds hsml = Float64(HSML[p])

        in_image = false

        @inbounds for dim = 1:3

            @inbounds in_image = check_in_image(pos[dim], hsml,
                                                minCoords[dim], maxCoords[dim])

            # exit the loop if the particle is not in the image frame
            if !in_image
                break
            end
        end

        # only calculate the properties if the particle is in the image
        if in_image

            particles_in_image += 1

            # save rest of variables
            @inbounds hsml_inv    = Float64(1.0/hsml)
            @inbounds bin_q       = Float64(Bin_Quant[p])
            @inbounds m           = Float64(M[p])
            @inbounds rho_inv     = Float64(1.0/ρ[p])

            pixmin = Vector{Int64}(undef,2)
            pixmax = Vector{Int64}(undef,2)

            @inbounds for dim = 1:2

                pixmin[dim] = find_min_pixel(pos[dim], hsml, minCoords[dim],
                                             param.pixelSideLength)

                pixmax[dim] = find_max_pixel(pos[dim], hsml, minCoords[dim],
                                             param.pixelSideLength, max_pixel[dim])

            end

            if conserve_quantities

                # calculate pixel area
                pixsize_inv = 1.0/param.pixelSideLength

                xp1 = (pos[1] - hsml) * pixsize_inv
                xp2 = (pos[1] + hsml) * pixsize_inv
                yp1 = (pos[2] - hsml) * pixsize_inv
                yp2 = (pos[2] + hsml) * pixsize_inv

                pix_area = (xp2 - xp1) * (yp2 - yp1) * param.pixelSideLength^2

                d3 = (m * rho_inv) / ( pix_area * param.pixelSideLength)

                # number of pixels over which the particle is distributed
                N_distr = Int64((pixmax[1] - pixmin[1] + 1 ) * (pixmax[2] - pixmin[2] + 1))

                # allocate arrays for weights
                kernel_tab = zeros(N_distr)
                d1_tab     = zeros(N_distr)
                d2_tab     = zeros(N_distr)
                dx_tab     = zeros(N_distr)

                # pixel count for attributing weights to pixels
                N_count = 1

                # first loop to calculate weights
                @inbounds for i = pixmin[2]:pixmax[2]
                    dy = param.y[i] - pos[2]
                    djmin = max(yp1, i - 1.0)
                    djmax = min(yp2, i)

                    @inbounds for j = pixmin[1]:pixmax[1]

                        dx = param.x[j] - pos[1]
                        dimin = max(xp1, j - 1.0)
                        dimax = min(xp2, j)

                        d1_tab[N_count] = dimax - dimin
                        d2_tab[N_count] = djmax - djmin

                        # compute distance to pixel center in units of hsml
                        dx_tab[N_count] = get_d_hsml_2D(dx, dy, hsml_inv)
                        # update pixel value
                        kernel_tab[N_count] = kernel_value_2D(kernel, dx_tab[N_count], hsml_inv)

                        N_count += 1

                    end # end x-loop
                end # end y-loop

                # reset the counter
                N_count = 1

            end # conserve_quantities


            # second loop to calculate value
            bin_prefac = bin_q * m * rho_inv

            @inbounds for i = pixmin[2]:pixmax[2]

                @inbounds for j = pixmin[1]:pixmax[1]


                    if !conserve_quantities

                        # calculate simple distance to pixel center
                        dx = param.x[j] - pos[1]
                        dy = param.y[i] - pos[2]
                        distance_hsml = get_d_hsml_2D(dx, dy, hsml_inv)

                        # update pixel value
                        image[j, i] += bin_prefac * kernel_value_2D(kernel, distance_hsml, hsml_inv)
                    else

                        # update pixel value with weights
                        image[j, i] += bin_prefac * d1_tab[N_count] * d2_tab[N_count] * d3 * kernel_tab[N_count]
                        N_count += 1
                    end

                end # end x-loop
            end # end y-loop

        end # end check if in image

        # update for ProgressMeter
        if show_progress
            #lock(P_lock)  # uncomment to make thread-safe if needed in the future
            idx += 1
            ProgressMeter.update!(P, idx)
            #unlock(P_lock)
        end
    end

    if show_progress
        @info "Mapped $particles_in_image particles."
    end

   return image
end


function sphMapping_3D(Pos, HSML, M, ρ, Bin_Quant;
                          param::mappingParameters, kernel::SPHKernel,
                          show_progress::Bool=true)

    N = length(M)  # number of particles

    image = zeros( length(param.x), length(param.y), length(param.z))

    minCoords = [param.x[1], param.y[1], param.z[1]]
    maxCoords = [param.x[end], param.y[end], param.z[end]]

    max_pixel = [length(param.x), length(param.y), length(param.z)]

    if show_progress
        P = Progress(N)
        idx = 0
        #P_lock = SpinLock()  # uncomment to make thread-safe if needed in the future
    end

    @inbounds for p = 1:N

        # save stuff from array to single variables
        pos      = Float64.(Pos[p,:])
        hsml     = Float64(HSML[p])

        in_image = false

        @inbounds for dim = 1:3

            in_image = check_in_image(pos[dim], hsml,
                                      minCoords[dim], maxCoords[dim])

            # exit the loop if the particle is not in the image frame
            if !in_image
                break
            end
        end

        if in_image

            # save rest of variables
            hsml_inv = Float64(1.0/hsml)
            bin_q    = Float64(Bin_Quant[p])
            m        = Float64(M[p])
            rho_inv  = Float64(1.0/ρ[p])

            pixmin = Vector{Int64}(undef,3)
            pixmax = Vector{Int64}(undef,3)

            @inbounds for dim = 1:3

                pixmin[dim] = find_min_pixel(pos[dim], hsml, minCoords[dim],
                                             param.pixelSideLength)

                pixmax[dim] = find_max_pixel(pos[dim], hsml, minCoords[dim],
                                             param.pixelSideLength, max_pixel[dim])

            end

            bin_prefac = bin_q * m * rho_inv

            @inbounds for k = pixmin[3]:pixmax[3]
                dz = param.z[k] - pos[3]

                @inbounds for j = pixmin[2]:pixmax[2]
                    dy = param.y[j] - pos[2]

                    @inbounds for i = pixmin[1]:pixmax[1]
                        dx = param.x[i] - pos[1]

                        # compute distance to pixel center in units of hsml
                        distance_hsml = get_d_hsml_3D(dx, dy, dz, hsml_inv)

                        # update pixel value
                        image[i, j, k] += bin_prefac * kernel_value_3D(kernel, distance_hsml, hsml_inv)

                    end # end x-loop
                end # end y-loop
            end # end z-loop

        end # end check if in image

        # update for ProgressMeter
        if show_progress
            #lock(P_lock)  # uncomment to make thread-safe if needed in the future
            idx += 1
            ProgressMeter.update!(P, idx)
            #unlock(P_lock)
        end
    end

   return image
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

                        norm += kernel_value_3D(kernel, distance/hsml, Float64.(hsml))

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

                        val[i,j,k] += bin_quant * mass / rho * kernel_value_3D(kernel, distance/hsml, Float64.(hsml))  * norm

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

include(joinpath(dirname(@__FILE__), "kernels.jl"))
include(joinpath(dirname(@__FILE__), "sph_types.jl"))

using Statistics
# function sphCenterMapping(part::Dict{Any,Any}, broadening::Float64,
#                           param::mappingParameters, kernel)
#
#     N = length(part["POS"][:,1])  # number of particles
#
#     global progress = 0
#     global done = 0
#
#     for p = 1:N
#
#         hsml::Float64 = part["HSML"][i] + broadening
#
#         pixmin = Vector{Int64}(undef,3)
#         pixmax = Vector{Int64}(undef,3)
#
#         minCoords::Vector{Float64}(undef,3) = [ pixelXCoordinates[1],
#                                                 pixelYCoordinates[1],
#                                                 pixelZCoordinates[1] ]
#         for i = 1:3
#
#             pixmin[i] = ( part["POS"][p,i] + hsml - minCoords[i] ) / param.pixelSideLength
#             if pixmin[i] < 0
#                 pixmin[i] = 0
#             end
#
#             pixmax[i] = ( ( part["POS"][p,i] + hsml - minCoords[i] ) / param.pixelSideLength + 1 )
#
#             max = ( i < 2 : param.pixelsPerDimension : zPixels)
#             if pixmax[i] > max
#                 pixmax[i] = max
#             end
#
#         end
#
#         if  pixmin[1] >= pixmax[1] ||
#             pixmin[2] >= pixmax[2] ||
#             pixmin[3] >= pixmax[3]
#
#             global progress += 1
#             #continue        # ???
#
#         end
#
#         for i = pixmin[1]:pixmax[1]
#
#             for j = pixmin[2]:pixmax[2]
#
#                 for k = pixmin[3]:pixmax[3]
#
#                     d = sqrt.( (pixelXCoordinates[i] - part["POS"][p,1])^2 +
#                                (pixelYCoordinates[i] - part["POS"][p,2])^2 +
#                                (pixelZCoordinates[i] - part["POS"][p,3])^2  )
#
#                     # woher kommt massLoc?
#                     sum = massLoc[p] / densityLoc[p] * kernel_value(kernel, d/hsml, hsml)
#
#                     if selection == 0 || selection[p]
#
#                         valuesLoc[k] += dataLoc[offsetInData + p * stepInData] * sum
#
#                     end
#
#                     normsLoc[k] += sum
#
#                 end
#
#             end
#
#         end
#
#     end
#
# end


function get_bounds(pos::Array{Float32,2}, min::Float64, max::Float64, axis::Int64,
                    range_arr::Vector{Int64})

    k = findall(pos[range_arr,axis] .>= min  )
    m = findall(pos[range_arr[k],axis] .<= max  )

    return range_arr[k[m]]

end
#
#
# function get_dist_to_center(cell::Vector{Float64}, pos::Array{Float32,1})
#
#     return sqrt.(abs(cell[1] - pos[1])^2 +
#                  abs(cell[2] - pos[2])^2 +
#                  abs(cell[3] - pos[3])^2 )
#
# end
#
# function get_kernel_values(pos::Array{Float32,2}, hsml::Array{Float32,1},
#                            cell::Vector{Float64}; kernel)
#
#     global val = 0.
#
#     for i = 1:length(pos[:,1])
#
#         d = get_dist_to_center(cell, pos[i,:])
#         global val += kernel_value(kernel, d, Float64(hsml[i]))
#
#     end
#
#     return val
#
# end
#
#
# function sphCenterMapping(part::Dict{Any,Any}, broadening::Float64,
#                           param::mappingParameters; kernel)
#
#
#     hsml = maximum(part["HSML"])
#
#     #println("hsml: ", hsml)
#
#     x = zeros(length(param.x)*length(param.y))
#     y = zeros(length(param.x)*length(param.y))
#     c = zeros(length(param.x)*length(param.y))
#
#     global total_cells = length(c)
#
#     range_arr = collect(1:length(part["POS"][:,1]))
#
#     #println("range: ", length(range_arr))
#
#     global cellnum = 1
#
#     # x-axis
#     for i ∈ 1:(length(param.x_max)-1)
#
#         range_arr_x = get_bounds(part["POS"],
#                            (param.x_max[i] - 2. * hsml),
#                            (param.x_max[i+1] + 2. * hsml),
#                            1, range_arr)
#
#         #println("x: ", length(range_arr_x))
#         # y-axis
#         for j ∈ 1:(length(param.y_max)-1)
#
#             range_arr_y = get_bounds(part["POS"],
#                                (param.y_max[j] - 2. * hsml),
#                                (param.y_max[j+1] + 2. * hsml),
#                                2, range_arr_x)
#
#             #println("y: ", length(range_arr_y))
#
#             if length(range_arr) > 0
#
#                 # z-axis
#                 for k ∈ 1:(length(param.z_max)-1)
#
#                     range_arr_z = get_bounds(part["POS"],
#                                        (param.z_max[k] - 2. * hsml),
#                                        (param.z_max[k+1] + 2. * hsml),
#                                        3, range_arr_y)
#
#                     if length(range_arr) > 0
#
#                         #println("z: ", length(range_arr_z))
#
#                         cell_center = [param.x[i], param.y[j], param.z[k]]
#
#                         c[cellnum] += get_kernel_values(part["POS"][range_arr_z,:],
#                                                         part["HSML"][range_arr_z],
#                                                         cell_center,
#                                                         kernel=kernel)
#                     end
#                 end
#             end
#
#             x[cellnum] = param.x[i]
#             y[cellnum] = param.y[j]
#             #println("Done with cell ", cellnum, " of ", total_cells)
#             #println("value c: ", c[cellnum])
#
#             global cellnum += 1
#
#         end
#
#         #println(100. * (cellnum/total_cells), "% done")
#     end
#
#
#     return x, y, c
#
#
# end


function sphAdaptiveMapping(Pos::Array{Float32,2}, HSML::Array{Float32,2},
                            M::Array{Float32,1}, ρ::Array{Float32,2},
                            Bin_Quant::Array{Float32,2},
                            broadening::Float64,
                            param::mappingParameters; kernel)

    #println("mapping ", length(total_range), " particles")

    N = length(M)
    pix_Vol = param.pixelSideLength^3
    boundary = false

    val = zeros(length(param.x), length(param.y),length(param.z))

    pixmin = Vector{Int64}(undef,3)
    pixmax = Vector{Int64}(undef,3)

    minCoords = [param.x[1], param.y[1], param.z[1]]

    for p = 1:N

        pos = Pos[p,:]
        hsml = HSML[p]
        mass = M[p]
        rho = ρ[p]
        bin_quant = Bin_Quant[p]

        for dim = 1:3

            pixmin[dim] = Int64( round((pos[dim] - hsml - minCoords[dim]) / param.pixelSideLength ))

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
            closestPixel = Vector{Int64}(undef,3)
            closestDistance = Inf

            for i = pixmin[1]:pixmax[1]
                for j = pixmin[2]:pixmax[2]
                    for k = pixmin[3]:pixmax[3]

                        distance = sqrt.( (param.x[i] - pos[1])^2 +
                                          (param.y[j] - pos[2])^2 +
                                          (param.z[k] - pos[3])^2 )

                        if distance < closestDistance
                            closestPixel = [i, j, k]            # ?
                            closestDistance = distance          # ?
                        end

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

                        distance = sqrt.( (param.x[i] - pos[1])^2 +
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

    cell_count = 1
    for i = 1:length(param.x)
        for j = 1:length(param.y)


            #c[i,j] = sum(val[i,j,:])
            c[j,i] = sum(val[i,j,:]) # transpose for imshow()

            # x[cell_count] = param.x[i]
            # y[cell_count] = param.y[j]
            # c[cell_count] = sum(val[i,j,:])
            #
            # cell_count += 1

        end
    end


    return c
#    return x, y, c

end

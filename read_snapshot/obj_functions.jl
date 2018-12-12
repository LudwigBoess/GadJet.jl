# ------------------------------------------------------------------------------
# --------------- Object functions ---------------------------------------------
# ------------------------------------------------------------------------------

function head_to_obj(filename)

    h = Header()

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

    h.npart = read!(f, Array{Int32,1}(undef,6))
    h.massarr = read!(f, Array{Float64,1}(undef,6))
    h.time = read(f, Float64)
    h.z = read(f, Float64)
    h.flag_sfr = read(f, Int32)
    h.flag_feedback = read(f, Int32)
    h.nall = read!(f, Array{UInt32,1}(undef,6))
    h.flag_cooling = read(f, Int32)
    h.num_files = read(f, Int32)
    h.boxsize = read(f, Float64)
    h.omega_0 = read(f, Float64)
    h.omega_l = read(f, Float64)
    h.h0 = read(f, Float64)
    h.flag_stellarage = read(f, Int32)
    h.flag_metals = read(f, Int32)
    h.npartTotalHighWord = read!(f, Array{UInt32,1}(undef,6))
    h.flag_entropy_instead_u = read(f, Int32)
    h.flag_doubleprecision = read(f, Int32)
    h.flag_ic_info = read(f, Int32)
    h.lpt_scalingfactor = read(f, Float32)

    close(f)

    return h

end


# function snap_to_obj(filename)
#
#     head = head_to_obj(filename)
#
#     if head.snap_format == 2
#
#        P = snap_2_p(filename, head)
#
#     else
#
#        P = snap_1_p(filename, head)
#
#     end
#
#     return P
#
# end


# function snap_2_p(filename, head)
#
#         f = open(filename)
#
#         seek(f, 296)
#
#         N = sum(head.npart)
#         skipsize = read(f, Int32, 1)[1]
#         bit_size = skipsize/(3*N)
#
#         if Int(bit_size) == 4
#             P = Array{part_single,1}(0)
#             # set up array of particles
#
#             for i = 1:6
#                  if head.npart[i] != Int32(0)
#                      for k = 1:head.npart[i]
#                          push!(P,part_single(i-1))
#                      end
#                  end
#             end
#
#         elseif Int(bit_size) == 8
#
#                 P = Array{part_double,1}(0)
#                 for i = 1:6
#                      if head.npart[i] != Int32(0)
#                          for k = 1:head.npart[i]
#                              push!(P,part_double(i-1))
#                          end
#                      end
#                 end
#         else
#             println("read error! neither 32 nor 64 bits data!")
#             return -1
#         end
#
#         # read positions
#         if Int(bit_size) == 4
#             for i = 1:sum(head.npart)
#                  P[i].pos = read(f, Float32, 3)
#             end
#         else
#             for i = 1:sum(head.npart)
#                  P[i].pos = read(f, Float64, 3)
#             end
#         end
#
#
#         p = position(f)
#
#         # skip identifiers
#         seek(f, p+24)
#
#         if Int(bit_size) == 4
#             dummy = read(f, Float32, (3,sum(head.npart)))
#         else
#             dummy = read(f, Float64, (3,sum(head.npart)))
#         end
#
#         # read Velocities
#         for i = 1:sum(head.npart)
#
#             P[i].vel = dummy[:,i]
#
#         end
#
#         p = position(f)
#
#         # skip identifiers
#         seek(f, p+24)
#
#         for i = 1:sum(head.npart)
#
#             P[i].id = read(f, UInt32)
#
#         end
#
#         p = position(f)
#
#         # skip identifiers
#         seek(f, p+24)
#
#         # read masses
#         for k = 1:6
#             if head.massarr[k] == Int32(0)
#                  for i = 1:head.npart[k]
#                      if bit_size == 4
#                         P[i].m = read(f, Float32)
#                      else
#                         P[i].m = read(f, Float64)
#                      end
#                  end
#             else
#                  for i = 1:head.npart[k]
#                     P[i].m = head.massarr[k]
#                  end
#             end
#         end
#
#
#         if head.npart[1] != Int32(0)
#
#             # skip identifiers
#             p = position(f)
#             seek(f, p+24)
#
#             # read U
#             for i = 1:head.npart[1]
#
#                if bit_size == 4
#                   P[i].U = read(f, Float32)
#                else
#                   P[i].U = read(f, Float64)
#                end
#
#            end
#
#             # skip identifiers
#             p = position(f)
#             seek(f, p+24)
#
#             # Read Density
#             for i = 1:head.npart[1]
#
#                if bit_size == 4
#                   P[i].rho = read(f, Float32)
#                else
#                   P[i].rho = read(f, Float64)
#                end
#
#            end
#
#             # skip identifiers
#             p = position(f)
#             seek(f, p+24)
#
#             # read smoothing length
#             for i = 1:head.npart[1]
#
#                if bit_size == 4
#                   P[i].sml = read(f, Float32)
#                else
#                   P[i].sml = read(f, Float64)
#                end
#
#            end
#
#
#         end
#
#         close(f)
#
#         return P
#
# end
#
# function snap_1_p(filename, head)
#
#     f = open(filename)
#
#     seek(f, 264)
#
#     N = sum(head.npart)
#     skipsize = read(f, Int32, 1)[1]
#     bit_size = skipsize/(3*N)
#
#     if Int(bit_size) == 4
#         P = Array{part_single,1}(0)
#         # set up array of particles
#
#         for i = 1:6
#              if head.npart[i] != Int32(0)
#                  for k = 1:head.npart[i]
#                      push!(P,part_single(i-1))
#                  end
#              end
#         end
#
#     elseif Int(bit_size) == 8
#
#             P = Array{part_double,1}(0)
#             for i = 1:6
#                  if head.npart[i] != Int32(0)
#                      for k = 1:head.npart[i]
#                          push!(P,part_double(i-1))
#                      end
#                  end
#             end
#     else
#         println("read error! neither 32 nor 64 bits data!")
#         return -1
#     end
#
#     # read positions
#     if Int(bit_size) == 4
#         for i = 1:sum(head.npart)
#              P[i].pos = read(f, Float32, 3)
#         end
#     else
#         for i = 1:sum(head.npart)
#              P[i].pos = read(f, Float64, 3)
#         end
#     end
#
#
#     p = position(f)
#
#     # skip identifiers
#     seek(f, p+8)
#
#     if Int(bit_size) == 4
#         dummy = read(f, Float32, (3,sum(head.npart)))
#     else
#         dummy = read(f, Float64, (3,sum(head.npart)))
#     end
#
#     # read Velocities
#     for i = 1:sum(head.npart)
#
#         P[i].vel = dummy[:,i]
#
#     end
#
#     p = position(f)
#
#     # skip identifiers
#     seek(f, p+8)
#
#     for i = 1:sum(head.npart)
#
#         P[i].id = read(f, UInt32)
#
#     end
#
#     p = position(f)
#
#     # skip identifiers
#     seek(f, p+8)
#
#     # read masses
#     for k = 1:6
#         if head.massarr[k] == Int32(0)
#              for i = 1:head.npart[k]
#                  if bit_size == 4
#                     P[i].m = read(f, Float32)
#                  else
#                     P[i].m = read(f, Float64)
#                  end
#              end
#         else
#              for i = 1:head.npart[k]
#                 P[i].m = head.massarr[k]
#              end
#         end
#     end
#
#
#     if head.npart[1] != Int32(0)
#
#         # skip identifiers
#         p = position(f)
#         seek(f, p+8)
#
#         # read U
#         for i = 1:head.npart[1]
#
#            if bit_size == 4
#               P[i].U = read(f, Float32)
#            else
#               P[i].U = read(f, Float64)
#            end
#
#        end
#
#         # skip identifiers
#         p = position(f)
#         seek(f, p+8)
#
#         # Read Density
#         for i = 1:head.npart[1]
#
#            if bit_size == 4
#               P[i].rho = read(f, Float32)
#            else
#               P[i].rho = read(f, Float64)
#            end
#
#        end
#
#         # skip identifiers
#         p = position(f)
#         seek(f, p+8)
#
#         # read smoothing length
#         for i = 1:head.npart[1]
#
#            if bit_size == 4
#               P[i].sml = read(f, Float32)
#            else
#               P[i].sml = read(f, Float64)
#            end
#
#        end
#
#
#     end
#
#     close(f)
#
#     return P
#
# end

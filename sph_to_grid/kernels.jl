mutable struct Cubic
    ngb::Int64
end

mutable struct Quintic
    ngb::Int64
end

mutable struct WendlandC4
    ngb::Int64
end

mutable struct WendlandC6
    ngb::Int64
end

function kernel_value(kernel::Cubic, u::Float64, h::Float64)

    norm = 8.0/π
    n = norm/h^3

    if u < 0.5
        return ( 1.0 + 6.0 * (u - 1.0) * u^2) * n
    elseif u < 1.0
        return ( 2.0 * (1.0 - u) * (1.0 - u) * (1.0 - u)) * n
    else
        return 0.
    end

end

function kernel_deriv(kernel::Cubic, u::Float64, h::Float64)

    norm = 8.0/π
    n = norm/h^4

    if u < 0.5
        return ( u * (18.0 * u - 12.0 )) * n
    elseif u < 1.0
        return ( -6.0 * (1.0 - u) * (1.0 - u) ) * n
    else
        return 0.
    end

end


function kernel_value(kernel::Quintic, u::Float64, h::Float64)

    norm = ( 2187.0 / ( 40. * π))
    n = norm/h^3

    if u < 1.0/3.0
        return ( ( 1.0 - u )^5 - 6.0 * ( 2.0/3.0 - u )^5  + 15.0 * ( 1.0/3.0 - u )^5 ) * n
    elseif u < 2.0/3.0
        return ( ( 1.0 - u )^5 - 6.0 * ( 2.0/3.0 - u )^5 ) * n
    elseif u < 1.0
        return ( ( 1.0 - u )^5 ) * n
    else
        return 0.
    end

end

function kernel_deriv(kernel::Quintic, u::Float64, h::Float64)

    norm = ( 2187.0 / ( 40. * π))
    n = norm/h^4

    if u < 1.0/3.0
        return ( -5.0 * ( 1.0 - u )^4 + 30.0 * ( 2.0/3.0 - u )^4  - 75.0 * ( 1.0/3.0 - u )^4 ) * n
    elseif u < 2.0/3.0
        return ( -5.0 * ( 1.0 - u )^4 + 30.0 * ( 2.0/3.0 - u )^4 - 75.0 ) * n
    elseif u < 1.0
        return ( -5.0 * ( 1.0 - u )^4 ) * n
    else
        return 0.
    end

end


function kernel_value(kernel::WendlandC4, u::Float64, h::Float64)

    norm = 495.0/(32. * π)
    n = norm/h^3

    if u < 1.0
        return ( ( 1. - u )^6 * ( 1.0 + 6.0 * u + 35.0/3.0 * u^2 ) ) * n
    else
        return 0.
    end

end

function kernel_deriv(kernel::WendlandC4, u::Float64, h::Float64)

    norm = 495.0/(32. * π)
    n = norm/h^4

    if u < 1.0
        return ( -288.0/3.0 * ( 1. - u )^5 * u^2 - 56.0/3.0 * u * ( 1. - u)^5 ) * n
    else
        return 0.
    end

end


function kernel_value(kernel::WendlandC6, u::Float64, h::Float64)

    norm = 1365.0/(64.0*π)
    n = norm/h^3

    if u < 1.0
        return ( (1.0 - u)^8 * ( 1.0 + 8. * u + 25. * u^2 + 32. * u^3 )) * n
    else
        return 0.
    end

end

function kernel_deriv(kernel::WendlandC6, u::Float64, h::Float64)

    norm = 1365.0/(64.0*π)
    n = norm/h^4

    if u < 1.0
        return ( -22. * (1.0 - u)^7 * u * ( 16. * u^2 + 7. * u + 1. )) * n
    else
        return 0.
    end

end



##################################
########## Sandbox ###############
##################################

# using PyPlot
# using BenchmarkTools
#
# function get_values()
#
#     u = collect(0.0:0.001:1.3)
#     h = 1.0
#
#     cu = Cubic(10)
#     qui = Quintic(10)
#     we4 = WendlandC4(10)
#     we6 = WendlandC6(10)
#
#     cu_val = zeros(length(u))
#     qui_val = zeros(length(u))
#     we4_val = zeros(length(u))
#     we6_val = zeros(length(u))
#
#     for i = 1:length(u)
#         cu_val[i] = kernel_value(cu, u[i], h)
#         qui_val[i] = kernel_value(qui, u[i], h)
#         we4_val[i] = kernel_value(we4, u[i], h)
#         we6_val[i] = kernel_value(we6, u[i], h)
#     end
#
#     return u, cu_val, qui_val, we4_val, we6_val
#
# end
#
#
# u, cu_val, qui_val, we4_val, we6_val = get_values()
#
#
# function plot_it(u, cu_val, qui_val, we4_val, we6_val)
#     fig = figure()
#     ax = gca()
#
#     ax[:set_ylim]([0,7.0])
#     ax[:set_xlim]([0.0,1.0])
#
#     plot(u, cu_val, label="Cubic")
#     plot(u, qui_val, label="Quintic")
#     plot(u, we4_val, label="Wendland C4")
#     plot(u, we6_val, label="Wendland C6")
#
#     xlabel("r/h")
#     ylabel("ρ")
#
#     title("SPH Kernels")
#
#     legend()
#
#     savefig("kernels.pdf")
#
#     close(fig)
#
# end
#
#
# plot_it(u, cu_val, qui_val, we4_val, we6_val)

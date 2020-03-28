"""
    calculate_center_of_mass(x, m)

Calculates the center of mass for given positions `x` and masses `m`.
"""
function calculate_center_of_mass(x, m)

    dt = typeof(m[1])
    com = Array{dt,1}(undef,3)

    com[1] = sum(m .* x[:,1])
    com[2] = sum(m .* x[:,2])
    com[3] = sum(m .* x[:,3])

    mtot = sum(m)

    com ./= mtot

    return com
end

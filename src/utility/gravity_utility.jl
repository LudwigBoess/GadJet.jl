"""
    calculate_center_of_mass(x, m)

Calculates the center of mass for given positions `x` and masses `m`.
"""
function calculate_center_of_mass(x, m)

    dt = typeof(m[1])
    com = zeros(dt, 3)

    for i = 1:length(m)
        com[1] += m[i] * x[i,1]
        com[2] += m[i] * x[i,2]
        com[3] += m[i] * x[i,3]
    end
    mtot = sum(m)

    com ./= mtot

    return com
end

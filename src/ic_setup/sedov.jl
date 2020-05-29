
function make_sedov(fi, file_out, kernel;
                    Ndim::Int64=3, code_units::Bool=false,
                    u_min::Float64=1.e-15, hdf5::Bool=false,
                    GU::GadgetPhysical=GadgetPhysical())

    h = head_to_obj(fi)

    pos_info  = Info_Line("POS", Float32, 3, [1,0,0,0,0,0])
    hsml_info = Info_Line("HSML", Float32, 1, [1,0,0,0,0,0])

    x0    = read_block_by_name(fi, "POS", info=pos_info, parttype=0)
    hsml0 = read_block_by_name(fi, "HSML", info=hsml_info, parttype=0)

    if Ndim == 3
        x0    .*= 6.0
        hsml0 .*= 6.0
    end

    #x, hsml = build_box(x0, 4, hsml0)
    x = x0
    #x = build_grid(50)
    hsml = hsml0

    if Ndim == 3
        x = x .- 3.0
    else
        x[:,1] .-= 0.5
        x[:,2] .-= 0.5
        x[:,3]  .= 0.0
    end

    r = sqrt.( x[:,1].^2 + x[:,2].^2 .+ x[:,3].^2)

    p_center = sortperm(r)[1]

    r_0 = sqrt.( (x[:,1] .- x[p_center,1]).^2 .+
                 (x[:,2] .- x[p_center,2]).^2 .+
                 (x[:,3] .- x[p_center,3]).^2   )

    hsml_center = Float64(hsml[p_center])

    i_energy = sortperm(r_0)[1:(kernel.n_neighbours+1)]

    weights = zeros(kernel.n_neighbours+1)

    if Ndim == 3
        for i = 1:kernel.n_neighbours+1
            u = r_0[i_energy[i]]/hsml_center
            weights[i] = kernel_value_3D(kernel, u, hsml_center)
        end
    else
        for i = 1:kernel.n_neighbours+1
            u = r_0[i_energy[i]]/hsml_center
            weights[i] = kernel_value_2D(kernel, u, hsml_center)
        end
    end

    # normalize weight
    weights ./= weights[1]
    weights ./= sum(weights)

    if code_units
        ρ = 1.0
        V = 1.0
    else
        #ρ = 1.24e7

        V = 6.0^3

        m_unit = 1.989e43
        l_unit = 3.085678e21
        mp = 1.6726219e-24

        ρ = (l_unit^3 * mp)/m_unit
    end

    M = ρ*V

    N = length(x[:,1])
    m = Float32.( (M/N) .* ones(N))

    println("Number of particles: ", N)
    println("Mass of single particle: ", m[1])

    u = zeros(Float32, N)

    if code_units
        E = 1.0/m[1]
    else
        E = (1.e51/GU.E_cgs)/m[1]
    end

    Tmax = E * GU.T_K
    Tmin = 10.0 # K

    u_min = Tmin / GU.T_K
    #E = Tmax / GU.T_K
    println("Maximum Temperature: $Tmax K")
    println("Minimum Temperature: $Tmin K")

    println("Minimum Energy: $u_min")

    println("Energy to distribute: $E")

    energies = Float32.(weights .* E)

    u[i_energy] .= energies

    u[u .< u_min] .= u_min

    println("Distributed energy: $(sum(u[i_energy]))")

    println("Error: $((sum(u[i_energy]) - E)/sum(u[i_energy]))")

    id = UInt32.(collect(0:N-1))

    v = [ Float32.(zeros(N)) Float32.(zeros(N)) Float32.(zeros(N)) ]

    B_seed = 1.e-12
    B = [ Float32.(B_seed .* ones(N)) Float32.(zeros(N)) Float32.(zeros(N))]

    if Ndim == 3
        x .+= 3.0
    else
        x[:,1] .+= 0.5
        x[:,2] .+= 0.5
        x[:,3]  .= 0.0
    end


    println("Setup done...")

    compare_h = head_to_obj(fi)

    head = compare_h
    head.time = 0.0
    head.npart[1] = N
    head.nall[1] = N
    head.massarr[1] = m[1]

    if Ndim == 3
        head.boxsize = 6.0
    else
        head.boxsize = 1.0
    end

    x = Float32.(x)
    hsml = Float32.(hsml)

    f = open(file_out, "w")
    write_header(f, head)
    write_block(f, x, "POS")
    write_block(f, v, "VEL")
    write_block(f, id, "ID")
    write_block(f, u, "U")
    write_block(f, hsml, "HSML")
    write_block(f, B, "BFLD")
    close(f)


    println()
end

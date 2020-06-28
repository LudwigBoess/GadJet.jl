using Printf
using LinearAlgebra: norm
using FFTW
using ProgressMeter

struct ShockParameters

    glass_file::String
    output_file::String
    n_blocks::Int64
    v::Array{Float64,2}
    B::Array{Float64,2}
    U::Vector{Float64}
    B0::Float64
    turb::Bool

    function ShockParameters(glass_file::String="", output_file::String="",
                             U::Vector{Float64}=zeros(2),
                             n_blocks::Int64=70,
                             v::Array{Float64,2}=zeros(2,3),
                             B::Array{Float64,2}=zeros(2,3),
                             B0::Float64=0.0, turb::Bool=false)

        new(glass_file, output_file, n_blocks,
            v, B, U, B0, turb)

    end
end

function P_to_U(P::Real, rho::Real)
    return P/((2.0/3.0) * rho)
end

function get_bfield_from_angle(theta::Float64;
                               reduction_scale::Float64=1.0e-10)

    shock_normal = [1.0, 0.0, 0.0]
    # convert to radians
    theta *= π/180.0

    bfld = zeros(3)

    bfld[1] = sin(90.0 * π/180.0 - theta)
    bfld[2] = sin(theta)

    bfld ./= sqrt(bfld[1]^2 + bfld[2]^2 + bfld[3]^2)

    #bfld += shock_normal

    bfld .*= reduction_scale

    return bfld

end

function getLargeBox(x, hsml=0)

    x = [x; x[:,1] .+ 1 x[:,2:3]]
    x = [x; x[:,1] x[:,2] .+ 1 x[:,3]]
    x = [x; x[:,1:2] x[:,3] .+ 1]

    x = x ./ 2

    if hsml != 0

        h = []
        for i = 1:8
            h = vcat(h, hsml)
        end

        return x, h

    else

        return x

    end
end

function buildTube(x0, n_blocks, hsml0=0, offset=0)

    n_part = length(x0[:,1])
    x = zeros(n_blocks*n_part, length(x0[1,:]))

    for i = 0:n_blocks-1
        x[i*n_part+1:(i+1)*n_part,:] = i .* [1.0 0.0 0.0] .+ x0 .+ [offset 0.0 0.0]
    end

    if hsml0 != 0
        hsml = zeros(n_blocks*n_part)
        for i = 0:n_blocks-1
            hsml[i*n_part+1:(i+1)*n_part] = hsml0
        end

        return Float32.(x), Float32.(hsml)
    end

    return Float32.(x)
end

function build_B_tube(B_in, n_blocks)

    n_part = length(B_in[:,1])
    B = zeros(n_blocks*n_part, length(B_in[1,:]))

    for i = 0:n_blocks-1
        B[i*n_part+1:(i+1)*n_part,:] .= B_in
    end

    return Float32.(B)
end

function findii(x, xtab)

    n = length(xtab)

    imin  = -1
    imax  = -1
    dmin1 = maximum(xtab) - minimum(xtab)
    dmin2 = maximum(xtab) - minimum(xtab)

    for i = 2:n
        d = x - xtab[i]
        if (0.0 < d < dmin1)
            imin = i
            dmin1 = d
        end

        if ( d < 0.0  && -d < dmin2 )
            imax = i
            dmin2 = -d
        end
    end
    return imin, imax
end

function setup_turb_B(pos, npart, B0)

    r = zeros(npart)

    for i = 1:length(r)
        r[i] = norm(pos[i,:])
    end

    k = 1.0./r

    # set grid
    # NFFT  = 117
    # NFFT2 = 58

    NFFT  = 233
    NFFT2 = 116

    minn = maximum(abs.(k)) / minimum(abs.(k))
    mink = minimum(abs.(k))

    k0 = maximum(abs.(k))

    α = 5.0/3.0

    A0 = k0^(-α) * B0^2

    α2 = 0.5*α

    nnn = (NFFT2+1)^3 * 4

    val = randn(nnn)
    ϕ2  = rand(nnn) .* 2.0π

    kx = zeros(Int64, NFFT)
    ky = zeros(Int64, NFFT)
    kz = zeros(Int64, NFFT)

    xp = zeros(Float64, NFFT)
    yp = zeros(Float64, NFFT)
    zp = zeros(Float64, NFFT)

    Bhx = zeros(NFFT, NFFT, NFFT)
    Bhy = zeros(NFFT, NFFT, NFFT)
    Bhz = zeros(NFFT, NFFT, NFFT)


    ii = 1

    @showprogress "Constructing grid..." for iz=0:NFFT2

        niz = (iz > 0) ? (NFFT-iz) : 0
        kz[iz+1]  =  iz
        kz[niz+1] = -iz

        for iy = 0:NFFT2

            niy = ( iy > 0 ) ? (NFFT - iy) : 0
            ky[iy+1]  =  iy
            ky[niy+1] = -iy

            for ix = 0:NFFT2

                nix = ( ix > 0 ) ? (NFFT - ix) : 0
                kx[ix+1]  =  ix
                kx[nix+1] = -ix

                k2 = ( kx[ix+1]^2 + ky[iy+1]^2 + kz[iz+1]^2 ) * mink^2
                PSP = ( k2 != 0.0 ) ? sqrt(A0 * k2^α2 ) :  0.0

                # first quadrant

                Bh = PSP * val[ii]
                if  kz[iz+1] == 0
                    ϕ1 = π/4.0
                else
                    ϕ1 = atan( - (kx[ix+1]*cos(ϕ2[ii])
                               +  ky[iy+1]*sin(ϕ2[ii]))
                               /  kz[iz+1] )
                end

                Bhx[ix+1, iy+1, iz+1]    = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[ix+1, iy+1, iz+1]    = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[ix+1, iy+1, iz+1]    = Bh*sin(ϕ1)

                Bhx[nix+1, niy+1, niz+1] = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[nix+1, niy+1, niz+1] = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[nix+1, niy+1, niz+1] = Bh*sin(ϕ1)


                # second quadrant
                ii += 1
                Bh = PSP * val[ii]
                if  kz[iz+1] == 0
                    ϕ1 = π/4.0
                else
                    ϕ1 = atan( - (kx[ix+1]*cos(ϕ2[ii])
                               +  ky[iy+1]*sin(ϕ2[ii]))
                               /  kz[iz+1] )
                end

                Bhx[nix+1, iy+1, iz+1]   = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[nix+1, iy+1, iz+1]   = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[nix+1, iy+1, iz+1]   = Bh*sin(ϕ1)

                Bhx[ix+1, niy+1, niz+1]  = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[ix+1, niy+1, niz+1]  = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[ix+1, niy+1, niz+1]  = Bh*sin(ϕ1)

                # third quadrant
                ii += 1
                Bh = PSP * val[ii]
                if  kz[iz+1] == 0
                    ϕ1 = π/4.0
                else
                    ϕ1 = atan( - (kx[ix+1]*cos(ϕ2[ii])
                               +  ky[iy+1]*sin(ϕ2[ii]))
                               /  kz[iz+1] )
                end

                Bhx[ix+1, niy+1, iz+1]   = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[ix+1, niy+1, iz+1]   = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[ix+1, niy+1, iz+1]   = Bh*sin(ϕ1)

                Bhx[nix+1, iy+1, niz+1]  = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[nix+1, iy+1, niz+1]  = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[nix+1, iy+1, niz+1]  = Bh*sin(ϕ1)

                # fourth quadrant
                ii += 1
                Bh = PSP * val[ii]
                if  kz[iz+1] == 0
                    ϕ1 = π/4.0
                else
                    ϕ1 = atan( - (kx[ix+1]*cos(ϕ2[ii])
                               +  ky[iy+1]*sin(ϕ2[ii]))
                               /  kz[iz+1] )
                end

                Bhx[nix+1, niy+1, iz+1]  = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[nix+1, niy+1, iz+1]  = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[nix+1, niy+1, iz+1]  = Bh*sin(ϕ1)

                Bhx[ix+1, iy+1, niz+1]   = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[ix+1, iy+1, niz+1]   = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[ix+1, iy+1, niz+1]   = Bh*sin(ϕ1)


                ii += 1
            end

        end

    end

    #Bhx
    Bhx_dummy = Bhx
    Bhx = Bhx_dummy

    for iz = 1:NFFT
        xp[iz] = 1.0/(kx[iz] * mink)
        yp[iz] = 1.0/(ky[iz] * mink)
        zp[iz] = 1.0/(kz[iz] * mink)
    end

    Bhx = fft(Bhx)
    Bfehler1 = sum(abs.(imag.(Bhx)))
    Bhx = real.(Bhx)
    Bfehler2 = sum(abs.(Bhx))

    Bhy = fft(Bhy)
    #Bfehler1 = sum(abs.(imag.(Bhx)))
    Bhy = real.(Bhy)
    #Bfehler2 = sum(abs.(Bhx))

    Bhz = fft(Bhz)
    #Bfehler1 = sum(abs.(imag.(Bhx)))
    Bhz = real.(Bhz)
    #Bfehler2 = sum(abs.(Bhx))

    anzz=NFFT^3
    bmean  = sqrt((sum(Bhx)/anzz)^2+(sum(Bhy)/anzz)^2+(sum(Bhz)/anzz)^2)
    bmean2 = sum(sqrt.(Bhx.^2+Bhy.^2+Bhz.^2))/anzz

    println("Gitter: <|B|>= $bmean2, |<B>|= $bmean")

    x = pos[:,1]
    y = pos[:,2]
    z = pos[:,3]

    bx = zeros(npart)
    by = zeros(npart)
    bz = zeros(npart)

    @showprogress "Interpolating B..." for i = 1:npart

        xvec = pos[i,:]

        ixmin, ixmax = findii(pos[i,1], xp)
        iymin, iymax = findii(pos[i,2], yp)
        izmin, izmax = findii(pos[i,3], zp)

        x1vec = zeros(3)
        x2vec = zeros(3)

        x1vec[1], x2vec[1]= xp[ixmin], xp[ixmax]
        x1vec[2], x2vec[2]= xp[iymin], xp[iymax]
        x1vec[3], x2vec[3]= xp[izmin], xp[izmax]

        dx = ( x[i] - xp[ixmin] ) / ( xp[ixmax] - xp[ixmin] )
        dy = ( y[i] - yp[ixmin] ) / ( yp[ixmax] - yp[ixmin] )
        dz = ( z[i] - zp[ixmin] ) / ( zp[ixmax] - zp[ixmin] )

        fx1y1z1 = Bhx[ixmin, iymin, izmin]
        fx1y1z2 = Bhx[ixmin, iymin, izmax]
        fx1y2z1 = Bhx[ixmin, iymax, izmin]
        fx1y2z2 = Bhx[ixmin, iymax, izmax]
        fx2y1z1 = Bhx[ixmax, iymin, izmin]
        fx2y1z2 = Bhx[ixmax, iymin, izmax]
        fx2y2z1 = Bhx[ixmax, iymax, izmin]
        fx2y2z2 = Bhx[ixmax, iymax, izmax]

        bx[i]= ( 1.0 - dx ) * ( 1.0 - dy )* ( 1.0 - dz ) * fx1y1z1 +
               ( 1.0 - dx ) * ( 1.0 - dy )*         dz   * fx1y1z2 +
               ( 1.0 - dx ) *         dy  * ( 1.0 - dz ) * fx1y2z1 +
               ( 1.0 - dx ) *         dy  *         dz   * fx1y2z2 +
                       dx   * ( 1.0 - dy) * ( 1.0 - dz)  * fx2y1z1 +
                       dx   * ( 1.0 - dy) *         dz   * fx2y1z2 +
                       dx   *         dy  * ( 1.0 - dz)  * fx2y2z1 +
                       dx   *         dy  *         dz   * fx2y2z2

        # y component
        fx1y1z1 = Bhy[ixmin, iymin, izmin]
        fx1y1z2 = Bhy[ixmin, iymin, izmax]
        fx1y2z1 = Bhy[ixmin, iymax, izmin]
        fx1y2z2 = Bhy[ixmin, iymax, izmax]
        fx2y1z1 = Bhy[ixmax, iymin, izmin]
        fx2y1z2 = Bhy[ixmax, iymin, izmax]
        fx2y2z1 = Bhy[ixmax, iymax, izmin]
        fx2y2z2 = Bhy[ixmax, iymax, izmax]

        by[i]= ( 1.0 - dx ) * ( 1.0 - dy )* ( 1.0 - dz ) * fx1y1z1 +
               ( 1.0 - dx ) * ( 1.0 - dy )*         dz   * fx1y1z2 +
               ( 1.0 - dx ) *         dy  * ( 1.0 - dz ) * fx1y2z1 +
               ( 1.0 - dx ) *         dy  *         dz   * fx1y2z2 +
                       dx   * ( 1.0 - dy) * ( 1.0 - dz)  * fx2y1z1 +
                       dx   * ( 1.0 - dy) *         dz   * fx2y1z2 +
                       dx   *         dy  * ( 1.0 - dz)  * fx2y2z1 +
                       dx   *         dy  *         dz   * fx2y2z2

        # z component
        fx1y1z1 = Bhz[ixmin, iymin, izmin]
        fx1y1z2 = Bhz[ixmin, iymin, izmax]
        fx1y2z1 = Bhz[ixmin, iymax, izmin]
        fx1y2z2 = Bhz[ixmin, iymax, izmax]
        fx2y1z1 = Bhz[ixmax, iymin, izmin]
        fx2y1z2 = Bhz[ixmax, iymin, izmax]
        fx2y2z1 = Bhz[ixmax, iymax, izmin]
        fx2y2z2 = Bhz[ixmax, iymax, izmax]

        bz[i]= ( 1.0 - dx ) * ( 1.0 - dy )* ( 1.0 - dz ) * fx1y1z1 +
               ( 1.0 - dx ) * ( 1.0 - dy )*         dz   * fx1y1z2 +
               ( 1.0 - dx ) *         dy  * ( 1.0 - dz ) * fx1y2z1 +
               ( 1.0 - dx ) *         dy  *         dz   * fx1y2z2 +
                       dx   * ( 1.0 - dy) * ( 1.0 - dz)  * fx2y1z1 +
                       dx   * ( 1.0 - dy) *         dz   * fx2y1z2 +
                       dx   *         dy  * ( 1.0 - dz)  * fx2y2z1 +
                       dx   *         dy  *         dz   * fx2y2z2
    end

    return [bx by bz]
end


function setup_shocktube(par::ShockParameters)

    println("reading glass file")
    pos_info = Info_Line("POS", Float32, 3, [1,0,0,0,0,0])
    hsml_info = Info_Line("HSML", Float32, 1, [1,0,0,0,0,0])

    h = head_to_obj(par.glass_file)

    pos = Float32(1.0/h.boxsize) .* read_block_by_name(par.glass_file, "POS", info=pos_info, parttype=0)
    hsml = Float32(1.0/h.boxsize) .* read_block_by_name(par.glass_file, "HSML", info=hsml_info, parttype=0)

    println("read!")

    println("setting up tube")
    x_large, hsml_l = getLargeBox(pos, hsml)

    m = 1.0/length(x_large[:,1])

    x_small = pos

    n_blocks = par.n_blocks

    # build tubes
    println("Building x tubes")
    x_left, hsml_left = buildTube(x_large, n_blocks, hsml_l)
    x_right, hsml_right = buildTube(x_small, n_blocks, hsml, n_blocks)

    x = [x_left; x_right]
    hsml_out = [hsml_left; hsml_right]
    N = length(x[:,1])

    left_part = findall(x[:,1] .<= Float64(par.n_blocks))
    right_part = findall(x[:,1] .> Float64(par.n_blocks))

    # set up random magnetic field
    if par.turb

        n_large = length(x_large[:,1])
        B_large = setup_turb_B(x_large, n_large, par.B0)

        n_small = length(x_small[:,1])
        B_small = setup_turb_B(x_small, n_small, par.B0)

        println("Building B tubes")
        B_left = build_B_tube(B_large, n_blocks)
        B_right = build_B_tube(B_small, n_blocks)
        B = [B_left; B_right]
        N_B = length(B[:,1])

        if N != N_B
            error("Error in tube building!\nN = $N\nN_B = $N_B")
        end

    else
         B = [ Float32.(zeros(N)) Float32.(zeros(N)) Float32.(zeros(N))]

        # Bx =  0.75  0.75
        # By = +1    -1
        # Bz =  0     0

        B[left_part,1]  .= Float32(par.B[1,1])
        B[right_part,1] .= Float32(par.B[2,1])
        B[left_part,2]  .= Float32(par.B[1,2])
        B[right_part,2] .= Float32(par.B[2,2])
        B[left_part,3]  .= Float32(par.B[1,3])
        B[right_part,3] .= Float32(par.B[2,3])
    end

    println("done")
    MASS = Float32.(m .* ones(N))

    println("Checking uniqueness of positions.")
    unique_check = (length(unique(x, dims=1)[:,1]) == N) ? true : false

    if unique_check == false
        return "ERROR! Overlapping particles!"
    end


    println(minimum(x[:,1]), " ", maximum(x[:,1]))

    println("Number of particles = ", N)

    println("Assigning shock parameters")


    head = head_to_obj(par.glass_file)
    head.boxsize = 100000.0
    head.npart[1] = N
    head.nall[1] = N
    head.massarr[1] = m

    ρ = Vector{Float32}(undef, N)
    ρ[left_part] .= Float32(1.0)
    ρ[right_part] .= Float32(0.125)

    vel = [ Float32.(zeros(N)) Float32.(zeros(N)) Float32.(zeros(N))]

    id = UInt32.(collect(0:N-1))

    U = Vector{Float32}(undef, N)

    U[left_part] .= Float32(par.U[1])
    U[right_part] .= Float32(par.U[2])

    println("done")

    println("writing ic file")

    f = open(par.output_file, "w")
    write_header(f, head)
    write_block(f, x, "POS")
    write_block(f, vel, "VEL")
    write_block(f, id, "ID")
    write_block(f, U, "U")
    write_block(f, hsml_out, "HSML")
    write_block(f, B, "BFLD")
    close(f)

end

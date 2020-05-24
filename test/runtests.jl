using GadJet, Test, DelimitedFiles, Unitful, UnitfulAstro

@testset "Read Snapshot" begin

    snap_file = joinpath(dirname(@__FILE__), "snap_050")

    @test_nowarn read_snap(snap_file)
    @test_nowarn read_header(snap_file)
    @test_nowarn read_info(snap_file)

    d = read_snap(snap_file, "POS", 0)

    ideal_file = joinpath(dirname(@__FILE__), "pos.dat")
    d_ideal = Float32.(readdlm(ideal_file))

    @test d == d_ideal
end

@testset "Unit Conversion" begin

    GU = GadgetPhysicalUnits()

    @test GU.t_s ≈ 3.085678e16u"s"
    @test GU.E_cgs ≈ 1.989e53u"erg"

    GU = GadgetPhysical()

    @test GU.t_s ≈ 3.085678e16
    @test GU.E_cgs ≈ 1.989e53
end

@testset "Riemann Sod-Shock" begin

    par = RiemannParameters(Ul=100.0, Mach=10.0, t=1.5)

    sol = solve([86.0], par)

    @test sol.rho3 ≈ 0.37667409437005994
    @test sol.rho4 ≈ 0.4854368932038835
    @test sol.P34  ≈ 13.097397114653086

end


@testset "Riemann CR-Shock" begin

    par = RiemannParameters(Pl=63.400, Pr=0.1, Mach=10.0, dsa_model=4, t=1.5)

    sol = solve([86.0], par)

    @test sol.rho3    ≈ 0.3699882303652922
    @test sol.rho4    ≈ 0.5925766595991485
    @test sol.P34_tot ≈ 12.089335761741005
    @test sol.P4_cr   ≈ 3.583255961157783

end

@testset "Synchrotron Kernel" begin

    @test synchrotron_kernel(  1.0) ≈ 0.6514228153553652
    @test synchrotron_kernel(  3.5) ≈ 0.08268719536694746
    @test synchrotron_kernel(  4.5) ≈ 0.03357441971502366
    @test synchrotron_kernel( 10.0) ≈ 0.00019223826430086885
    @test synchrotron_kernel( 50.0) ≈ 1.734785203976593e-21
    @test synchrotron_kernel(100.0) ≈ 4.697593665922202e-43

end

@testset "SPH Kernels" begin

    k = Cubic()
    d = kernel_value_2D(k, 0.5, 0.5)
    @test d ≈ 0.15915494309189535
    d = kernel_value_3D(k, 0.5, 0.5)
    @test d ≈ 0.07957747154594767

    k = Quintic()
    d = kernel_value_2D(k, 0.5, 0.5)
    @test d ≈ 0.07767855829318415
    d = kernel_value_3D(k, 0.5, 0.5)
    @test d ≈ 0.06630419797168217

    k = WendlandC4()
    d = kernel_value_2D(k, 0.5, 0.5)
    @test d ≈ 0.30960610023345264
    d = kernel_value_3D(k, 0.5, 0.5)
    @test d ≈ 0.26606774238812336

    k = WendlandC6()
    d = kernel_value_2D(k, 0.5, 0.5)
    @test d ≈ 0.052822211162893276
    d = kernel_value_3D(k, 0.5, 0.5)
    @test d ≈ 0.05055250677698771

end


@testset "glimpse" begin

    snap_file = joinpath(dirname(@__FILE__), "snap_050")
    d = glimpse(snap_file, "RHO", [3.0, 3.0, 3.0], 6.0, 6.0, 6.0,
                parallel=false, verbose=false)

    ideal_file = joinpath(dirname(@__FILE__), "image.dat")
    d_ideal = readdlm(ideal_file)

    @test d[  1,  1] ≈ d_ideal[1, 1]
    @test d[ 30, 32] ≈ d_ideal[30, 32]
    @test d[117, 92] ≈ d_ideal[117, 92]

end


@testset "COM" begin

    m = ones(5)

    x = [  3.2  4.1 -2.4;
           2.9 -2.3  1.7;
          -2.3  4.2  3.2;
          -1.1 -2.4 -1.1;
           2.8 -3.2 -3.2 ]

    result = calculate_center_of_mass(x, m)

    @test result[1] ≈  1.1
    @test result[2] ≈  0.07999999999999999
    @test result[3] ≈ -0.36000000000000004

end


@testset "Spectral CRs types" begin

    d = CRShockData(1)
    @test d.dt[1] ==  0.0

    d = CRMomentumDistributionConfig()
    @test d.mp == 1.0

    d = CRMomentumDistribution(1)
    @test d.CRp_dis[1] == 0.0

end

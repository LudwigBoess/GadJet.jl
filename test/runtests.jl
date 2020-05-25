using GadJet, Test, DelimitedFiles, Unitful, UnitfulAstro

@testset "GadJet" begin

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

    @testset "Write Snapshot" begin

        # read in reference file
        ref_file = joinpath(dirname(@__FILE__), "snap_050")
        head = head_to_obj(ref_file)
        x = read_snap(ref_file, "POS", 0)

        # specify output file for testing
        output_file = joinpath(dirname(@__FILE__), "write_test.dat")
        f = open(output_file, "w")

        @test_nowarn write_header(f, head)
        @test_nowarn write_block(f, x, "POS")
        @test_nowarn write_block(f, x, "", snap_format=1)

        @test_warn "Please specify blockname!" write_block(f, x, "")

        close(f)
    end

    @testset "Unit Conversion" begin

        GU = GadgetPhysicalUnits()

        @test GU.t_s ≈ 3.085678e16u"s"
        @test GU.E_cgs ≈ 1.989e53u"erg"

        d = strip_unit(1.0u"g")
        @test d == 1.0

        GU = GadgetPhysical()

        @test GU.t_s ≈ 3.085678e16
        @test GU.E_cgs ≈ 1.989e53


    end

    @testset "Riemann Parameters" begin

        @test_warn "Both Ul and Pl are zero!" SodCRParameters_noCRs()
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

    @testset "DSA models" begin
        @test_warn "Invalid DSA model selection!" RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=5, t = 1.5)

        # KR07
        par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=0, t = 1.5)
        @test par.acc_function(5.0) ≈ 0.25185919999999995
        # KR13
        par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=1, t = 1.5)
        @test par.acc_function(5.0) ≈ 0.09999999999999998
        # Ryu+19
        par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=2, t = 1.5)
        @test par.acc_function(5.0) ≈ 0.017286554080677037
        # CS14
        par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=3, t = 1.5)
        @test par.acc_function(5.0) ≈ 0.04999999999999999
        # Pfrommer+16
        par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=4, t = 1.5)
        @test par.acc_function(5.0) == 0.5
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
        # < 0.5
        d = kernel_value_2D(k, 0.4, 0.5)
        @test d ≈ 0.26992678348385446
        d = kernel_value_3D(k, 0.4, 0.5)
        @test d ≈ 0.13496339174192723
        # < 1.0
        d = kernel_value_2D(k, 0.5, 0.5)
        @test d ≈ 0.15915494309189535
        d = kernel_value_3D(k, 0.5, 0.5)
        @test d ≈ 0.07957747154594767
        # > 1.0
        d = kernel_value_2D(k, 1.5, 0.5)
        @test d == 0.0
        d = kernel_value_3D(k, 1.5, 0.5)
        @test d == 0.0


        k = Quintic()
        # < 1/3
        d = kernel_value_2D(k, 0.3, 0.5)
        @test d ≈ 0.32700352517410614
        d = kernel_value_3D(k, 0.3, 0.5)
        @test d ≈ 0.2791208661307549
        # 2/3
        d = kernel_value_2D(k, 0.5, 0.5)
        @test d ≈ 0.07767855829318415
        d = kernel_value_3D(k, 0.5, 0.5)
        @test d ≈ 0.06630419797168217
        # < 1.0
        d = kernel_value_2D(k, 0.8, 0.5)
        @test d ≈ 0.0008155658657050455
        d = kernel_value_3D(k, 0.8, 0.5)
        @test d ≈ 0.0006961437210839494
        # > 1.0
        d = kernel_value_2D(k, 1.5, 0.5)
        @test d == 0.0
        d = kernel_value_3D(k, 1.5, 0.5)
        @test d == 0.0



        k = WendlandC4()
        # < 1.0
        d = kernel_value_2D(k, 0.5, 0.5)
        @test d ≈ 0.30960610023345264
        d = kernel_value_3D(k, 0.5, 0.5)
        @test d ≈ 0.26606774238812336
        # > 1.0
        d = kernel_value_2D(k, 1.5, 0.5)
        @test d == 0.0
        d = kernel_value_3D(k, 1.5, 0.5)
        @test d == 0.0




        k = WendlandC6()
        # < 1.0
        d = kernel_value_2D(k, 0.5, 0.5)
        @test d ≈ 0.052822211162893276
        d = kernel_value_3D(k, 0.5, 0.5)
        @test d ≈ 0.05055250677698771
        # > 1.0
        d = kernel_value_2D(k, 1.5, 0.5)
        @test d == 0.0
        d = kernel_value_3D(k, 1.5, 0.5)
        @test d == 0.0

    end


    @testset "SPH mappingParameters" begin

        @test_warn "Giving a center position requires extent in x, y and z direction." mappingParameters()

        @test_warn "Please specify pixelSideLenght or number of pixels!" mappingParameters(center=[0.0, 0.0, 0.0]
                                                                                           xlim = [-1.0, 1.0],
                                                                                           ylim = [-1.0, 1.0],
                                                                                           zlim = [-1.0, 1.0])

        @test_nowarn mappingParameters(center=[0.0, 0.0, 0.0]
                                       xlim = [-1.0, 1.0],
                                       ylim = [-1.0, 1.0],
                                       zlim = [-1.0, 1.0],
                                       pixelSideLength=0.2)
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

    @testset "Smac utility" begin

        @test_nowarn write_smac1_par("", 0, "", "", "", "",
                                     0, 0, 0, 4, 3,
                                     20.0, 10.0, 1, 1,
                                     24, 1.0, 1.e6, 10,
                                     1, 0.0, 0.0, 0.0)

        @test_warn "Read error: Incorrect image format!" read_smac1_binary_image(joinpath(dirname(@__FILE__), "image.dat"))

        filename = joinpath(dirname(@__FILE__), "GadJet_test.140.a.z.pix"))
        info = read_smac1_binary_info(filename)

        @test info.snap == 140

        image = read_smac1_binary_image(filename)
        @test length(image[:,1]) == 128
        @test image[1,1] ≈ 0.000120693

        @test_nowarn write_smac2_par(1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                     1.0, 1.0, 1024,
                                     "", "", "")
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

end

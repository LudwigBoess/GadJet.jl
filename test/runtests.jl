using GadJet, Test, DelimitedFiles

@testset "Read Snapshot" begin

    @info "Testing snapshot reading."

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

    @info "Testing unit conversion."

    GU = GadgetPhysicalUnits()

    @test GU.t_s ≈ 3.085678e16
    @test GU.E_cgs ≈ 1.989e53
end


@testset "Riemann Sod-Shock" begin

    @info "Testing Riemann solver for hydro Sod-Shock."

    par = RiemannParameters(Ul=100.0, Mach=10.0, t=1.5)

    sol = solve([86.0], par)

    @test sol.rho3 ≈ 0.37667409437005994
    @test sol.rho4 ≈ 0.4854368932038835
    @test sol.P34  ≈ 13.097397114653086

end


@testset "Riemann CR-Shock" begin

    @info "Testing Riemann solver for CR-Shock with acceleration."

    par = RiemannParameters(Pl=63.400, Pr=0.1, Mach=10.0, dsa_model=4, t=1.5)

    sol = solve([86.0], par)

    @test sol.rho3    ≈ 0.3699882303652922
    @test sol.rho4    ≈ 0.5925766595991485
    @test sol.P34_tot ≈ 12.089335761741005
    @test sol.P4_cr   ≈ 3.583255961157783

end

# @testset "glimpse" begin
#
#     @info "Testing glimpse function"
#
#     snap_file = joinpath(dirname(@__FILE__), "snap_050")
#     d = glimpse(snap_file, "RHO", [3.0, 3.0, 3.0], 6.0, 6.0, 1.0)
#
#     ideal_file = joinpath(dirname(@__FILE__), "image.dat")
#     d_ideal = readdlm(ideal_file)
#
#     @test d ≈ d_ideal
#
# end

include("wave_tests.jl")
using Test
tol = 10.0^(-5)
@testset "2D" begin
    @testset "Interior problems" begin
        @test closed2D() < tol
        @test cup2D() < tol
        @test open2D() < tol
    end
    @testset "Exterior problems" begin
        @test cyl_scat() < tol
        @test const2Dcylinder() < tol
        @test Hconst2Dcylinder() < tol
        @test nurbs2Dcylinder() < tol
    end
end

@testset "3D" begin
    @testset "Interior problems" begin
        @test closed3D() < tol
        @test cup3D() < tol
        @test open3D() < tol
    end
    @testset "Exterior problems" begin
        @test cylinder3D() < tol
        @test sphere() < tol
        @test cyl_scat3D() < tol
        @test sph_scat() < tol
    end
end

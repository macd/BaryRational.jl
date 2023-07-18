using LinearAlgebra
using BaryRational
using Test
using Random

include("test_aaa.jl")
include("test_fh.jl")
include("test_bary.jl")

@testset "BaryRational.jl" begin
    @testset "aaa_rational" begin
        @test test_aaa_gamma_poles()
        @test test_aaa_spiral()
        @test test_aaa_exp()
        @test test_aaa_length()
        @test test_aaa_tan()
        @test test_aaa_gamma()
        @test test_aaa_ranlog()
        @test test_aaa_infinp()
        @test test_aaa_naninp()
        @test test_aaa_residuals()
        #@test test_aaa_case2()
        #@test test_aaa_scale_invar()
        @test test_aaa_full_svd()
        @test test_aaa_deriv()
        @test test_aaa_deriv2()
        @test test_aaa_maxiters()
    end
    @testset "FH_rational_interpolation" begin
        @test test_fh_runge()
        @test test_abs_x()
        @test test_fh_deriv()
    end
    @testset "Barycentric interpolation on the Chebyshev points" begin
        @test test_bary_airy()
        @test test_bary_airy(BigFloat, tol=BigFloat(1//10^39))
    end
end

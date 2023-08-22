using LinearAlgebra
using BaryRational
using Test
using Random

include("test_aaa.jl")
include("test_fh.jl")
include("test_bary.jl")
include("test_deriv.jl")

@testset "BaryRational.jl" begin
    @testset "aaa_rational" begin
        @test test_aaa_gamma_poles()
        @test test_aaa_spiral()
        # These can take some time
        #@test test_big_aaa_spiral(verbose=true, clean=0) == false
        #@test test_big_aaa_spiral(verbose=true, clean=1)
        #@test test_big_aaa_spiral(verbose=true, clean=2)
        #@test test_big_aaa_spiral(verbose=true, clean=3)
        @test test_aaa_exp()
        @test test_aaa_length()
        @test test_aaa_tan()
        @test test_aaa_gamma()
        @test test_aaa_ranlog()
        @test test_aaa_infinp()
        @test test_aaa_naninp()
        @test test_aaa_residuals()
        #@test test_aaa_case2()
        @test test_aaa_scale_invar()
        @test test_aaa_full_svd()       # this has one doublet
        @test test_aaa_maxiters()       # 2 doublets
        @test test_aaa_truncation()
        @test test_aaa_complex()
        @test test_aaa_float32()
    end
    @testset "FH_rational_interpolation" begin
        @test test_fh_runge()
        @test test_abs_x()
        @test test_fh_complex()
    end
    @testset "Barycentric interpolation on the Chebyshev points" begin
        @test test_bary_airy()
        @test test_bary_airy(BigFloat, tol=BigFloat(1//10^39))
    end
    @testset "Test of Derivatives" begin
        @test test_aaa_deriv()
        @test test_aaa_deriv_pole()
        @test test_fh_deriv()
        @test test_aaa_airy_prime()
        @test test_2nd_derivative()
        @test test_runge_derivs()
        @test test_aaa_complex_deriv()
    end
end

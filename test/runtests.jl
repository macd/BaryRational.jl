using LinearAlgebra
using BaryRational
using Test
using Random

include("test_aaa.jl")
include("test_fh.jl")

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
    end
    @testset "FH_rational_interpolation" begin
        @test test_fh_runge()
        @test test_abs_x()
    end
end

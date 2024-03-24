# Borrowed from https://github.com/complexvariables/RationalFunctionApproximation.jl

pts = 10 .^ range(-15, 0, 500)
pts = [-reverse(pts); 0; pts]

approx(f; kw...) = aaa(f; kw...)

pass(f, r, z; kw...) = isapprox(f.(z), r(z), norm=u->norm(u,Inf); kw...)

@testset "Continuum basic functions" begin
    f = x -> abs(x + 0.5 + 0.01im); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> sin(1/(1.05-x)); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> exp(-1/(x^2)); @test pass(f, approx(f), pts, rtol=4e-13)
    f = x -> exp(-100x^2); @test pass(f, approx(f), pts, rtol=2e-13)
    f = x -> exp(-10/(1.2-x)); @test pass(f, approx(f), pts, rtol=1e-12)
    f = x -> 1/(1+exp(100*(x+.5))); @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> sin(100x) * exp(-10x^2); @test pass(f, approx(f), pts, atol=1e-11)
    f = x -> abs(x);  @test pass(f, approx(f), pts, atol=1e-8)
    f = x -> abs(x - 0.95);  @test pass(f, approx(f), pts, atol=1e-6)
end

@testset "Low-accuracy" begin
    f = x -> exp(3x);
    @test !pass(f, approx(f, tol=1e-4), pts, atol=1e-8)
end

@testset "Continuum Poles, zeros, residues" begin
    f = z -> (z+1) * (z+2) / ((z+3) * (z+4))
    r = approx(f)
    pol, res, zer = prz(r)
    @test isapprox(sum(pol+zer), -10, atol=1e-12)

    f = x -> 2/(3 + x) + 5/(x - 2im)
    r = approx(f)
    pol, res, zer = prz(r)
    @test isapprox(prod(res), 10, atol=1e-8)

    f = x -> sinpi(10x)
    r = approx(f)
    pol, res, zer = prz(r)
    @test isapprox(sort(abs.(zer))[19], 0.9, atol=1e-12)

    f = z -> (z - (3 + 3im))/(z + 2)
    r = approx(f)
    pol, res, zer = prz(r)
    @test isapprox(pol[1]*zer[1], -6-6im, atol=1e-12)
end

@testset "Continuum Vertical scaling" begin
    f = x -> 1e100*sin(x);
    @test pass(f, approx(f), pts, rtol=2e-13)
    f = x -> 1e-100*cos(x)
    @test pass(f, approx(f), pts, rtol=2e-13)
    #f = x -> 1e100*sin(x); @test pass(f, approx(f, =5, lawson=20), pts, rtol=1e-6)
    #f = x -> 1e-100*cos(x); @test pass(f, approx(f, =5, lawson=20), pts, rtol=1e-6)
end

@testset "Continuum Polynomials and reciprocals" begin
    args = Dict(:mmax=>150, :tol=>1e-13)
    f = x -> 0
    #@test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> x
    @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1im*x
    @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> x + x^2
    @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> x + x^3
    @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1/(1.1 + x)
    @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1/(1 + 1im*x)
    @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1/(3 + x + x^2)
    @test pass(f, approx(f; args...), pts, atol=2e-13)
    f = x -> 1/(1.01 + x^3)
    @test pass(f, approx(f; args...), pts, atol=2e-13)
end

# Note that mmax needs to be max_degree+1
@testset "Continuum Specified" begin
    f = x -> 0
    #@test pass(f, approx(f, mmax=1), pts, atol=2e-13)
    f = x -> x
    @test pass(f, approx(f, mmax=2), pts, atol=2e-13)
    f = x -> 1im*x
    @test pass(f, approx(f, mmax=4), pts, atol=2e-13)
    f = x -> x+x^2
    @test pass(f, approx(f, mmax=3), pts, atol=2e-13)
    f = x -> x+x^3
    @test pass(f, approx(f, mmax=4), pts, atol=2e-13)
    f = x -> 1/(1.1+x)
    @test pass(f, approx(f, mmax=4), pts, atol=2e-13)
    f = x -> 1/(1+1im*x)
    @test pass(f, approx(f, mmax=4), pts, atol=2e-13)
    f = x -> 1/(3+x+x^2)
    @test pass(f, approx(f, mmax=3), pts, atol=2e-13)
    f = x -> 1/(1.01+x^3)
    @test pass(f, approx(f, mmax=4), pts, atol=2e-13)
    f = x -> tanh(100x)
    #@test pass(f, approx(f), pts, atol=2e-13)
    f = x -> tanh(100*(x-.2))
    @test pass(f, approx(f), pts, atol=2e-13)
    @test pass(f, approx(f), pts, atol=2e-10)  # FIXME loss of accuracy
    f = x -> exp(x)
    @test pass(f, approx(f, tol=1e-13), pts, atol=2e-13)
    f = x -> cis(x)
    @test pass(f, approx(f), pts, atol=2e-13)
    f = x -> cis(x)
    r = approx(f);
    @test sum(@. abs(r(pts))-1) / length(pts) < 2e-13
    f = x -> exp(exp(x))/(x - 0.2im)
    r = approx(f)
    poles, res, zeros = prz(r)
    @test minimum(abs.(poles .- .2im)) < 1e-12
end

@testset "Continuum BigFloat" begin
    funcs = (cos, sin, exp)
    big_tol = BigFloat(1//BigInt(10)^70)
    for f in funcs
        @test pass(f, aaa(f, BigFloat), BigFloat.(pts), atol=big_tol)
    end
end

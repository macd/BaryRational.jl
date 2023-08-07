using SpecialFunctions
using ForwardDiff

# sort by increasing absolute value of the real part. Hack to get poles ordered
# for testing
csort(x) = sort(x, lt = (x,y) -> abs(real(x)) < abs(real(y)))

function test_aaa_spiral(T=Float64, verbose=false)
    zz = exp.(range(-0.5, complex(0.5, 15*pi), length=1000))

    yy = tan.(pi*zz/2)
    f = aaa(zz, yy, clean=false, verbose=verbose)
    pol, res, zer = prz(f)
    spol = csort(pol)
    verbose && foreach(println, real(spol))

    p1 = isapprox(abs(spol[1]),  1.0, atol=1e-14)
    p2 = isapprox(abs(spol[2]),  1.0, atol=1e-14)
    p3 = isapprox(abs(spol[3]),  3.0, atol=1e-6)
    p4 = isapprox(abs(spol[4]),  3.0, atol=1e-6)
    return all((p1, p2, p3, p4))
end

# So this turns out to be a good test for the various cleaning proceedures.
# It clean is turned off, we get spurious crap and the test will fail. For
# the other three cleaning routines, this did find some bugs. We do need to go
# to mmax=40 to see this. Can take a little bit of time.
function test_big_aaa_spiral(;verbose=true, clean=1)
    T = BigFloat
    zz = exp.(range(T(-1//2), complex(T(1//2), T(15//1)*pi), length=1000))
    yy = tan.(T(pi) *zz / T(2))
    f = aaa(zz, yy, mmax=40, tol=BigFloat(-1), # negative value forces mmax iters
            clean=clean, verbose=verbose)
    pol, res, zer = prz(f)
    spol = csort(pol)

    verbose && foreach(println, real(spol))

    p = falses(12)
    # This is actually amazing to me.
    ptols = [T(1//BigInt(10)^75), T(1//BigInt(10)^75),  # poles 1,-1
             T(1//BigInt(10)^55), T(1//BigInt(10)^55),  # poles 3,-3
             T(1//BigInt(10)^38), T(1//BigInt(10)^38),  # poles 5,-5
             T(1//BigInt(10)^28), T(1//BigInt(10)^28),  # poles 7,-7
             T(1//BigInt(10)^20), T(1//BigInt(10)^20),  # poles 9,-9
             T(1//BigInt(10)^15), T(1//BigInt(10)^15)]  # poles 11,-11
    for i in 1:2:12
        p[i]   = isapprox(abs(spol[i]),    T(i), atol=ptols[i])
        p[i+1] = isapprox(abs(spol[i+1]),  T(i), atol=ptols[i+1])
    end
    verbose && foreach(println, zip(1:length(p), p))
    return all(p)
end


# we test the values of the first 4 located poles to the accuracy noted in
# the AAA paper.  These values can make for brittle testing.
function test_aaa_gamma_poles()
    xx = range(-1.5, 1.5, length=100)
    yy = gamma.(xx)
    f = aaa(xx, yy)
    pol, res, zer = prz(f)
    spol = csort(pol)
    p1 = isapprox(spol[1],  0.0, atol=1e-15)
    p2 = isapprox(spol[2], -1.0, atol=1e-14)
    p3 = isapprox(spol[3], -2.0, atol=2e-7)
    p4 = isapprox(spol[4], -3.0, atol=3e-3)
    return p1 & p2 & p3 & p4
end


function test_aaa_exp(;tol=1e4*eps(1.0))
    Z = range(-1, 1, length=1000)
    F = exp.(Z)
    r = aaa(Z, F)
    t1 = norm(F - r(Z), Inf) < tol
    t2 = isnan(r(NaN))                        # check that r(NaN) = NaN
    t3 = !isinf(r(Inf))                       # r(inf) = sum(w.*f)/sum(w)
    return all((t1, t2, t3))
end

function test_aaa_length(tol=1e-8)
    Z = range(-1, 1, length=1000)
    F = exp.(Z)
    r = aaa(Z, F)
    m1 = length(r.x)
    r = aaa(Z, F, mmax=m1-1)
    t1 = length(r.x) == m1 - 1

    r = aaa(Z, F, tol=1e-3)
    t2 = length(r.x) < m1
    return all((t1, t2))
end

function test_aaa_tan(tol=1e-8)
    Z = range(-1, 1, length=1000)
    F = z -> tan(pi*z)
    r = aaa(Z, F)
    pol, res, zer  = prz(r)
    t1 = norm(F.(Z) - r(Z), Inf) < 10*tol
    t2 = minimum(abs.(zer)) < tol
    t3 = minimum(abs.(pol .- 0.5)) < tol
    t4 = minimum(abs.(res)) > 1e-13     # Test for spurious poles.
    return all((t1, t2, t3, t4))
end

function test_aaa_gamma(tol=1e-8)
    r = aaa([-0.9:0.05:1.0;], gamma)
    return abs(r(1.5) - gamma(1.5)) < 1e-3
end

function test_aaa_ranlog(tol=1e-8)
    Random.seed!(1137)
    Z = randn(10000) + 3im*randn(10000)
    f = z -> log(5 - z) / (1 + z^2)
    r = aaa(Z, f.(Z))
    return abs(r(0) - f(0)) < tol
end

function test_aaa_infinp(tol=1e-8)
    Z = complex.(range(-1.0, 1.0, length=101), zeros(101))
    r = aaa(Z, gamma.(Z))
    return abs(r(0.63) - gamma(0.63)) < 1e-3
end

function test_aaa_naninp(tol=1e-8)
    X = range(0, 20, length=100)
    F = sin.(X) ./ X
    r = aaa(X, F)
    return abs(r(2) - sin(2)/2) < 1e-3
end

function test_aaa_residuals(tol=1e-8)
    X = range(-1.337, 2, length=537)
    r = aaa(X, exp.(X) ./ X)
    pol, res, zer = prz(r)
    ii = findfirst(abs.(pol) .< 1e-8)
    t1 = abs(res[ii] - 1.0) < 1e-10

    r = aaa(X, (1+1im) .* gamma.(X))
    pol, res, zer = prz(r)
    ii = findfirst(abs.(pol .- (-1)) .< 1e-8)
    t2 = abs(res[ii] + (1+1im)) < 1e-10

    return all((t1, t2))
end


# The following test does not pass
function test_aaa_case2(tol=1e-8)
    # Case |Z| = 2: needs special treatment.
    Z = [0., 1]
    F = [1., 2]
    r = aaa(Z, F)
    t1 = norm(F .- r(Z), Inf) < tol
    t2 = r(Inf) == -Inf

    return all((t1, t2))
end

function test_aaa_scale_invar()
    # Check for exact scale-invariance
    Z = range(BigFloat(3//10), BigFloat(15//10), length=100)
    F = exp.(Z) ./ (BigFloat(1)+BigFloat(1)*im)
    # r1 = aaa(Z, F, clean=0)
    # r2 = aaa(Z, BigFloat(2)^311*F, clean=0)
    # r3 = aaa(Z, BigFloat(2)^-311*F, clean=0)
    r1 = aaa(Z, F)
    r2 = aaa(Z, BigFloat(2)^311*F)
    r3 = aaa(Z, BigFloat(2)^-311*F)
    t1 = r1(0.2im) == BigFloat(2)^-311*r2(0.2im)
    t2 = r1(1.4) == BigFloat(2)^311*r3(1.4)
    return all((t1, t2))
end


# Now a couple of tests of 'normal' functions on regular grids.

function do_aaa_func(func=abs, inc=0.01)
    x = [-1.0:inc:1.0;]
    y = func.(x)
    f = aaa(x, y)
    xt = [(-1.0 + inc):inc:(1.0 - inc);]
    yt = f.(xt)
    yexact = func.(xt)
    err = maximum(abs.(yt - yexact))
end

function test_aaa_abs_x(tol=1e-10)
    err = do_aaa_func()
    if err > tol
        println("AAA test failed  error: ", err, "   tol: ", tol)
        return false
    end
    return true
end

function runge_a(x)
    return 1.0 /(1 + x^2)
end


function test_aaa_runge(tol=1e-10)
    err = do_aaa_func(runge_a)
    if err > tol
        println("AAA test failed  error: ", err, "   tol: ", tol)
        return false
    end
    return true
end


# This setup requires that full=true for the svd so this tests that logic
# We also run both cleanup proceedures here.
function test_aaa_full_svd()
    x = [-1.0:0.2:1.0;]
    y = cos.(x)
    g  = aaa(x, y, clean=1)
    g2 = aaa(x, y, clean=2)
    return  isfinite(g(g.x[1])) &&  g.errvec[end] < 1e-14
end


# We approximate a function with a branch cut which requires a lot of
# poles/iterations since AAA clusters poles near branch points. See PR #8
function test_aaa_maxiters()
    f(a) = (q = sqrt(Complex(a^2 - 1)); (abs(q-a) <= 1 ? 1 : -1) * 2pi * inv(q))
    f(x, η) = f(cos(x) + im*η)
    ntrain = 10^4
    x = 2pi*range(0, step=1//ntrain, length=ntrain)
    z = cis.(x)
    eta = 1e-3
    fz = -im .* f.(x, eta) ./ z
    try
        aaa(z, fz, clean=1, verbose=false)
        aaa(z, fz, clean=2, verbose=false)
        return true
    catch e
        return !(e isa MethodError)
    end
end

function test_aaa_truncation()
    xbig = BigFloat.([-1//1:1//100:1//1;])
    fbig = sin.(xbig);
    # The min error appears at m=25, so this will test the truncation
    sf = aaa(xbig, fbig, mmax=30, clean=0, tol=BigFloat(1/10^40));
    # This also tests at support points
    xtest = BigFloat.([-1//1:1//1000:1//1;])
    error = norm(sin.(xtest) - sf.(xtest), Inf)
    return error < BigFloat(1//10^30)
end

function test_aaa_complex()
    x = [-1.0:0.1:1.0;]
    z = complex.(x,x)
    f = sin.(z)
    g = aaa(z, f)
    xx = [-1.0:0.001:1.0;]
    zz = complex.(xx, xx)
    return norm(sin.(zz) - g(zz), Inf) < 1e-12
end

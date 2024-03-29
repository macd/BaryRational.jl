using ForwardDiff
using ForwardDiff: derivative
using Random

function test_aaa_deriv()
    f(x) = sin(x) ^ 2 + cos(10x)
    x = [-1.0:0.01:1.0;]
    y = f.(x)
    g = aaa(x, y)
    Random.seed!(1137)
    xx = -1.0 .+ 2rand(100)
    dya = deriv.(g, xx)
    dyf = ForwardDiff.derivative.(f, xx)
    return  norm(dya - dyf, Inf) < 1e-10
end


# This one is challenging because we will have a downward spike at zero (and
# hence a discontinuity in the derivative) and we will have accuracy problems
# there. But elsewhere it is OK.
#
# TODO: rework this to use exponential clustering near zero using logrange for
#       the discrete AAA version. Also try the continuum AAA version
function test_aaa_deriv_pole()
    h(x) = cos(x) ^ 3 + sqrt(abs(sin(x)))
    x = [-1.0:0.01:1.0;]
    y = h.(x)
    g = aaa(x, y)
    xx = x .+ 0.005
    dya = deriv.(g, xx)
    dyf = ForwardDiff.derivative.(h, xx)
    idx = 94:106  # here be errors
    dya[idx] .= dyf[idx] .= 0.0
    return  norm(dya - dyf, Inf) < 1e-9
end


function test_fh_deriv()
    fh(x) = sin(x) ^ 2 + cos(10x)
    x = [-1.0:0.01:1.0;]
    y = fh.(x)
    fha = FHInterp(x, y, order=8, grid=true)
    xx = [-1.0:0.001:1.0;]
    dya = deriv.(fha, xx)
    dyf = ForwardDiff.derivative.(fh, xx)
    return  norm(dya - dyf, Inf) < 1e-8
end
    
# Test on the interval [-10.0, 0.0] where airyai is oscillatory
# and yet too small for asymptotic formulas to work.
function test_aaa_airy_prime(T=Float64; tol=T(1//10^10))
    xx = T(5) * (chebpts(256, T) .- T(1)) # move chebpts to [-10, 0.0] interval
    f = airyai.(xx)
    a = aaa(xx, f)

    # random test points
    Random.seed!(1137)
    xr = rand(T(-10):T(1//100):T(0), 1000)
    ybp = deriv.(a, xr)

    yap = airyaiprime.(xr)
    err1 = norm(ybp - yap, Inf)

    ap(x)  = derivative(airyai, x)
    ap2(x) = derivative(ap, x)

    err2 = norm(ap2.(xr) - deriv.(a, xr, m=2), Inf)
    return err1 < tol && err2 < tol
end

function test_2nd_derivative()
    f   = x ->  cos(10x)*exp(-x)
    df  = x -> -(10.0*sin(10x) + cos(10x))*exp(-x)
    df2 = x ->  (-99*cos(10x) + 20.0 * sin(10x)) * exp(-x)
    # The support points will be chosen from xx
    xx = [-1.0:0.01:1.0;]
    yy = f.(xx);
    g = aaa(xx, yy, do_sort=true)

    # These will naturally hit the support points too.
    xxx = [-1.0:0.001:1.0;]  

    err1 = norm(df.(xxx) - deriv.(g, xxx), Inf)
    err2 = norm(df2.(xxx) - deriv.(g, xxx, m=2), Inf)
    return err1 < 1e-9 && err2 < 1e-7
end

function test_runge_derivs(tol=1e-10)
    f(x) = 1.0 / (1.0 + x^2)
    xx = [-5.0:0.05:5.0;]
    yy = f.(xx)
    raaa = aaa(xx, yy, do_sort=true)

    yfd(x)  = derivative(raaa, x)
    yfd2(x) = derivative(yfd, x)
    yfd3(x) = derivative(yfd2, x)
    
    Random.seed!(1137)
    xr = 5.0 * (2rand(100) .- 1.0)
    yr = f.(xr)

    err1 = norm(yfd.(xr) - deriv.(raaa, xr), Inf)
    err2 = norm(yfd2.(xr) - deriv.(raaa, xr, m=2), Inf)
    err3 = norm(yfd3.(xr) - deriv.(raaa, xr, m=3), Inf)
    
    return err1 < tol && err2 < tol && err3 < 1e-5
end

# Not very accurate, but it seems to work.
function test_aaa_complex_deriv()
    x = [-1.0:0.01:1.0;]
    z = complex.(x, x)
    f = sin.(z)
    g = aaa(z, f)
    xx = [-1.0:0.001:1.0;]
    zz = complex.(xx, xx)
    return norm(cos.(zz) - deriv.(g, zz), Inf) < 1e-4
end

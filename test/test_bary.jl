# Test the simple barycentric interpolation on Chebyshev points.

# Test on the interval [-10.0, 0.0] where airyai is oscillatory
# and yet too small for asymptotic formulas to work.
function test_bary_airy(T=Float64; tol=T(1//10^13))
    xx = T(5) * (chebpts(256, T) .- T(1)) # move chebpts to [-10, 0.0] interval
    f = airyai.(xx)
    xr = rand(T(-10):T(1//100):T(0), 1000)
    yb = bary.(xr, (f,), (xx,))
    ya = airyai.(xr)
    err = abs.(yb - ya)
    return maximum(err) < tol
end

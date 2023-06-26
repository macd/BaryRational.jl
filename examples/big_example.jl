using ArbNumerics
using BaryRational

# The input grid to aaa. The support points will be picked from here
# in an adaptive fashion.
xx = [BigFloat(-10):BigFloat(1//100):BigFloat(0);];
yy = airyai.(xx);

# For BigFloat we need to set clean=false, I know from previous runs that
# mmax=53 gives the best result.
fa = aaa(xx, yy, clean=false, mmax=53, verbose=true, tol=BigFloat(1//10^40))

# Generate a 1000 random test points. Notice we take a much smaller step size 
# here to avoid landing on one of the support points
xrand = rand(BigFloat(-10):BigFloat(1//10000):BigFloat(0), 1000);

err = abs.(airyai.(xrand) - fa.(xrand));

println("Max error: ", maximum(err))

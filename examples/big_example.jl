using ArbNumerics
using BaryRational

# The input grid to aaa. The support points will be picked from here
# in an adaptive fashion.
xx = [BigFloat(-10):BigFloat(1//100):BigFloat(0);];
yy = airyai.(xx);

# For BigFloat we need to set clean=0, I know from previous runs that
# mmax=53 gives the best result but we let it run to 60 to test truncation code.
# Note that a negative tol forces mmax iterations
fa = aaa(xx, yy, clean=0, mmax=60, verbose=true, tol=BigFloat(-1))

# Generate a 1000 random test points. Notice we take a much smaller step size 
# here to (somewhat) avoid landing on one of the support points
xrand = rand(BigFloat(-10):BigFloat(1//10000):BigFloat(0), 1000);

err = norm(airyai.(xrand) - fa(xrand), Inf);

println("Max error: ", err)

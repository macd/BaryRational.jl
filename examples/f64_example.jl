using BaryRational
using SpecialFunctions

# The input grid to aaa. The support points will be picked from here
# in an adaptive fashion.
xx = [-10.0:0.01:0.0;];
yy = airyai.(xx);

fa = aaa(xx, yy, clean=false, verbose=true)

# Generate a 1000 random test points. Notice we take a much smaller step size 
# here to (somewhat) avoid landing on one of the support points
xrand = rand(-10.0:0.0001:0.0, 1000);

err = norm(airyai.(xrand) - fa.(xrand), Inf);

println("Max error: ", err)

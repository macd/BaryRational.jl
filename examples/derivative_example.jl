# Example of using deriv and how it gives correct results at a support point.
using PyPlot
using BaryRational
using ForwardDiff
using ForwardDiff: derivative

# The function and the true derivatives
f   = x ->  cos(10x)*exp(-x);
df  = x -> -(10.0*sin(10x) + cos(10x))*exp(-x);
df2 = x ->  (-99*cos(10x) + 20.0 * sin(10x)) * exp(-x);

# The support points will be chosen from xx
xx = [-1.0:0.01:1.0;]
yy = f.(xx);

g = aaa(xx, yy, do_sort=true, clean=1)

@show sw_deriv = deriv(g, 0.5)
@show df(0.5)
@show sw_deriv2 = deriv(g, 0.5, m=2)
@show df2(0.5)

gp1(x) = derivative(g,   x);
gp2(x) = derivative(gp1, x);

xxx = [-1.0:0.001:1.0;];

df_gold  = df.(xxx);
df2_gold = df2.(xxx);

err1fd = df_gold - gp1.(xxx);
err1sw = df_gold - deriv.(g, xxx);

@show norm(err1fd, Inf)
@show norm(err1sw, Inf)

fd2 = gp2.(xxx);

# Replace the values at support points in the ForwardDifference calculation
# with the true values of the 2nd derivative. Note that the points close
# to the support points will have increased error in the ForwardDifference
# calculation because the incorrect zeros in the 1st derivative propagate
# out from there (I think) and we make no attempt to correct for that.
II = findall(x -> x in g.x, xxx)
fd2[II] .= df2.(xxx)[II]

err2fd = df2_gold - fd2;
err2sw = df2_gold - deriv.(g, xxx, m=2);

@show norm(err2fd, Inf)
@show norm(err2sw, Inf)

plot(xxx, err2fd, color="red",   label="ForwardDifference 2nd derivative")
plot(xxx, err2sw, color="green", label="Scheider & Werner 2nd derivative")
xlabel("location")
ylabel("Error")
title("Error in 2nd derivative of rational approx to cos(10x)*exp(-x)")
legend()

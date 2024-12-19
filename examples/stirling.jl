using SpecialFunctions
using BaryRational

# This is just a fun example of "reproducing" Stirling results in
# interpolating the factorial. In Stirling's book "Methodus Differentialis"
# (see Ian Tweddle's translation) he interpolates the factorial function
# (in Proposition 21, example 2) to discover that (0.5)! = sqrt(pi) / 2.
# Holy heck, how did he do that? First, he uses his interpolation formula
# and goes to the tenth order (only even orders). Second, he actually
# interpolates the log of the factorial (who better than Stirling?) to find
# log((10.5)!)) and then exponentiates that and back tracks to 0.5.  Very cool.
# Here we use AAA instead of his interpolation formula.
function main()
    xx = [5:16;]
    yy = loggamma.(xx .+ 1)
    f = aaa(xx, yy)
    
    # t will be an approximation to (10.5)!
    t = exp(f(10.5))

    # Now we have (0.5)! = (10.5)! / 10.5 / 9.5 / 8.5 / ... / 1.5
    for i in 10:-1:1  
        t /= i + 0.5
    end

    sqpi2 = sqrt(pi) / 2
    println(t, "  ", sqpi2, "  ", abs(t - sqpi2))

    # But what if we ask AAA for log(0.5) directly? Well, we only get
    # a couple of digits. This is probably because we are extrapolating.
    s = exp(f(0.5))
    println(s, "  ", sqpi2, "  ", abs(s - sqpi2))

    # Let's give AAA a bit more to work with. So, yes, we do better,
    # about 6 digits this time
    xx = [0:16;]
    yy = loggamma.(xx .+ 1)
    f = aaa(xx, yy)
    s = exp(f(0.5))
    println(s, "  ", sqpi2, "  ", abs(s - sqpi2))

    # It makes one wonder if there is a general method to scale to a region
    # where AAA (or any approximation method) is more accurate, solve in
    # that region and then reverse the scaling to achieve a better final
    # accuracy.
end

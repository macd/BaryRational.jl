# NB: "This is not an item of software—it is a proof of concept."
# So take warning. This is brittle software. Change the sample points or
# even change the ordering of the sample points and results can
# drastically change. If you are trying to use BigFloat, you will be in
# for an adventure.

using BaryRational
using Base.MathConstants
using PyPlot
using SpecialFunctions
using ArbNumerics

#
# NOTE: For BigFloat computation we need Complex{BigFloat} for GenericShur.
#       You'd think we'd be always able to just return Complex{T} _but_
#       that gives worse results. See comment above.
#
function logspace(n, m, l; T=Float64)
    T <: Float64 && return T(10) .^ T.(range(n//1, m//1, length=l))
    Complex{T}(10) .^ Complex{T}.(range(n//1, m//1, length=l))
end

# Well, should we get rid of poles below or above the real axis? Depends on
# who you ask. If asking Born & Wolf, we need analytic and regular in the
# lower half plane.  If asking Morse & Feshbach, we need analytic and 
# regular in the upper half plane. We go for lower half here because that
# matches the sign in the Weideman paper, but it looks like Costa and Trefethen
# go with the upper half plane. Remarkably, except for a sign change, the
# results are very consistent when choosing either the upper or lower half plane.
# Also, not thrilled with the api for v,f that requires a vector.
function hilbert(u::Function; n=10, l=300, T=Float64, tol=T(1//10^13), clean=0,
                 verbose=false, mmax=100)
    X = logspace(-n, n, l; T=T)
    # If you did 'X = [-reverse(X); X]' instead, you get much worse results on
    # the Weideman testcases.
    X = [X; -X]     
    g = aaa(X, u.(X), clean=clean, tol=tol, verbose=verbose, mmax=mmax)
    pol, _, _ = prz(g)
    deleteat!(pol, imag(pol) .< T(0))
    verbose && println("Final pol length: ", length(pol))
    d = minimum(abs.(X .- transpose(pol)), dims=1)
    A = d ./ (X .- transpose(pol))
    A = [real(A) -imag(A)]
    c = reshape(A \ u.(X), :, 2) * Complex{T}.([1, 1im])
    # f is the complex analytic signal
    f = x -> reshape((d ./ (x[:] .- transpose(pol))) * c, size(x))
    # and its real part is the Hilbert transform of the real signal.
    v = x -> imag(f(x))                                  
    return v, f
end

# The functions and their Hilbert transforms are from "Computing the Hilbert 
# Transform on the Real Line" by J.A.C. Weideman
wfuncs = [x -> 1 / (1 + x^2),
          x -> 1 / (1 + x^4),
          x -> sin(x) / (1 + x^2),
          x -> sin(x) / (1 + x^4),
          x -> exp(-x^2),
          x -> sech(x),
          x -> exp(-abs(x))
]

# So dawson does not have a BigFloat or a Complex{BigFloat} version and
# digamma does not have Complex{BigFloat} version (it has a BigFloat though)
whilb = [
    y -> -y /(1 + y^2),
    y -> -y*(1+y^2)/(sqrt(2)*(1+y^4)),
    y -> (cos(y) - 1/e) / (1 + y^2),
    y -> (cos(y) - exp(-1/sqrt(2))*cos(1/sqrt(2)) - exp(-1/sqrt(2))*sin(1/sqrt(2))*y^2) / (1 + y^4),
    y -> -2 * dawson(y) / sqrt(pi),
    y -> tanh(y) + (im / pi) * (digamma(0.25 + y*im/(2pi)) - digamma(0.25 - y*im/(2pi))),
    y -> - (sign(y)/pi)*(exp(abs(y))*expint(abs(y)) + exp(-abs(y))*expinti(abs(y)))
]

function main(;T=Float64, clean=0, mmax=100, verbose=false, tol=T(1//10^13))
    for i in eachindex(wfuncs)
        v, f = hilbert(wfuncs[i]; T=T, verbose=verbose, mmax=mmax, tol=tol, clean=clean)
        err = abs(v([T(2)])[1] - whilb[i].(T(2)))
        println(v([T(2)])[1], "   ", whilb[i](T(2)), "   ", err)
    end
    v, f = hilbert(wfuncs[end], clean=clean)
    xx = range(T(-5), T(5), length=1000)
    yy = (v(xx) - whilb[end].(xx))
    plot(xx, yy)
    nothing
end


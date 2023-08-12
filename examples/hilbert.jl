# NB: "This is not an item of software—it is a proof of concept."
# So take warning. This is brittle software. Change the sample points or
# even change the ordering of the sample points and results can
# drastically change. If you are trying to use BigFloat, you will be in
# for an adventure (although some of it actually works OK)

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

hvec(x) = typeof(x) <: AbstractArray ? vec(x) : [x]

# Well, should we get rid of poles below or above the real axis? Depends on
# who you ask. If asking Born & Wolf, we need holomorphic in the
# lower half plane.  If asking Morse & Feshbach, we need holomorphic 
# in the upper half plane. We go for lower half here because that
# matches the sign in the Weideman paper, but it looks like Costa and Trefethen
# go with the upper half plane. Except for a sign change, the
# results are very consistent when choosing either the upper or lower half plane.
# Also, we follow the DSP.jl convention and only return the complex analytic signal,
# but here we return a function rather than the values of the Hilbert transform on
# the sample grid
function hilbert(X::AbstractVector{T}, Y::AbstractVector{T};
                 tol=T(1//10^13), clean=0, verbose=false, mmax=100) where {T}
    g = aaa(X, Y, clean=clean, tol=tol, verbose=verbose, mmax=mmax)
    pol, _, _ = prz(g)

    # Remove poles in the lower half plane
    deleteat!(pol, imag(pol) .< T(0))
    verbose && println("Final pol length: ", length(pol))

    d = minimum(abs.(X .- transpose(pol)), dims=1)
    A = d ./ (X .- transpose(pol))
    A = [real(A) -imag(A)]
    c = reshape(A \ Y, :, 2) * Complex{T}.([1, 1im])

    # f is a function that returns the complex analytic signal. It takes either
    # a scalar or vector, but (at present) always returns a vector.
    f = x -> reshape((d ./ (hvec(x) .- transpose(pol))) * c, size(x))

    return f
end


# If we are given a function to transform, follow C&T:
# "the sampling grid has been taken as 300 points exponentially spaced from 10^−10
# to 10^10 and their negatives, so 600 points all together."  Here we make these
# choices congifurable but defaulted to the values in the paper.
function hilbert(u::Function; n=10, l=300, T=Float64, tol=T(1//10^13), clean=0,
                 verbose=false, mmax=100)
    X = logspace(-n, n, l; T=T)
    X = [X; -X]   # This is what Costa and Trefethen do.
    #X = [-X; X]  # This is OK as well.
    #X = [-reverse(X); X]  # This is not. Much poorer results on Weideman examples
    return hilbert(X, u.(X), tol=tol, clean=clean, verbose=verbose, mmax=mmax)
end

# The functions and their Hilbert transforms are from "Computing the Hilbert 
# Transform on the Real Line" by J.A.C. Weideman, Mathematics of Computation (1995)
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
    f = nothing
    for i in eachindex(wfuncs)
        f = hilbert(wfuncs[i]; T=T, verbose=verbose, mmax=mmax, tol=tol, clean=clean)
        ht(x) = imag(f(x))[1]
        # These are the values and errors at x = 2.0 to replicate Costa & Trefethen
        err = abs(ht(T(2)) - whilb[i].(T(2)))
        println(ht(T(2)), "   ", whilb[i](T(2)), "   ", err)
    end

    # Plot the error of the Hilbert transform of the last test function
    xx = [-T(5):T(1//100):T(5);]
    yy = (imag(f(xx)) - whilb[end].(xx))
    plot(xx, yy);
end


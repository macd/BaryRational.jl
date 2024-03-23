# Is z within an isapprox of any x[i] ? (Only for use when we hit a NaN)
# Assumes that x is sorted from lowest to highest, which is kind of an
# imposition, but if not, then we are into linear search
function nearby(z::T, x::AbstractVector{T}) where {T}
    ord = T <: Complex ? cplxord : identity
    top = lastindex(x)
    bot = firstindex(x)
    mid = div(top, 2)
    while bot < mid < top
        if ord(z) > ord(x[mid])
            bot = mid
        elseif ord(z) < ord(x[mid])
            top = mid
        else
            return mid
        end
        mid = bot + div(top-bot, 2)
    end
    isapprox(x[bot], z) && return bot
    isapprox(x[top], z) && return top
    println(z, "  ", x[bot], "  ", x[top])
    error(z, " not found in x")
end


# For general location points with precalculated weights
""" 
    bary(z, f, x, w)  

evaluate f(z)

# Arguments
- z::Z                 : the point at which to evaluate f
- f::AbstractVector{F} : vector of function values at x
- x::AbstractVector{X} : vector of eval locations of f (sorted)
- w::AbstractVector{W} : weights for locations x
"""
function bary(z::Z, f::AbstractVector{F}, x::AbstractVector{X},
              w::AbstractVector{W}) where {Z, F, X, W}
    T = promote_type(Z, F, X, W)
    # assert(length(f) == length(x) == length(w))
    num = zero(T)
    den = zero(T)
    @inbounds for j in eachindex(f)
        t = w[j] / (z - x[j])
        num += t * f[j]
        den += t
    end
    fz = num / den
    fz = isfinite(fz) ? fz : f[nearby(z, x)]
end


function bary(z::T, a::U) where {T, U <: BRInterp}
    bary(z, a.f, a.x, a.w)
end


# When we don't get weigts we assume that x are at the Chebyshev points
# over the interval [x[begin], x[end]] and that w, the weights, are implied.
function bary(z, f::AbstractVector{T}, x) where {T}
    n = length(f)
    num = den = zero(z)
    t = T(1) / (T(2)*(z - x[1]))
    num += t * f[1]
    den += t
    sgn = T(-1)
    @inbounds for j = 2:n-1
        t = sgn / (z - x[j])
        num += t * f[j]
        den += t
        sgn = -sgn
    end
    t = sgn / (T(2) * (z - x[n]))
    num += t * f[n]
    den += t
    fz = num / den
    fz = isfinite(fz) ? fz : f[nearby(z, x)]
end


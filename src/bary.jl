# Is z within an isapprox of any x[i] ? 
# Assumes that x is sorted from lowest to highest
function nearby(z::T, x::Vector{T}) where {T}
    top = length(x)
    bot = 1
    mid = div(top, 2)
    while bot < mid < top
        if z > x[mid]
            bot = mid
        elseif z < x[mid]
            top = mid
        else
            return mid
        end
        mid = bot + div(top-bot, 2)
    end
    isapprox(x[bot], z) && return bot
    isapprox(x[top], z) && return top
    return -1   # should error
end


# For general location points with precalculated weights
""" 
    bary(z, f, x, w)  

evaluate f(z)

# Arguments
- `z::Float64`:         the point at which to evaluate f
- `f::Vector{Float64}`: vector of function values at x
- `x:Vector{Float64}`:  vector of eval locations of f (sorted)
- `w:Vector{Float64}`:  weights for locations x
"""
function bary(z::T, f::Vector{T}, x::Vector{T}, w::Vector{T}) where {T}
    # assert(length(f) == length(x) == length(w))
    num = zero(T)
    den = zero(T)
    @inbounds for j = 1:length(f)
        t = w[j] / (z - x[j])
        num += t * f[j]
        den += t
    end
    fz = num / den
    fz = isnan(fz) ? f[nearby(z, x)] : fz
end


function bary(z, f)
    n = length(f)
    x = chebpts(n)
    num = den = zero(z)
    t = 1 / 2(z - x[1])
    num += t * f[1]
    den += t
    sgn = -1
    @inbounds for j = 2:n-1
        t = sgn / (z - x[j])
        num += t*f[j]
        den += t
        sgn = -sgn
    end
    t = sgn / 2(z - x[n])
    num += t * f[n]
    den += t
    fz = num / den
    fz = isnan(fz) ? f[nearby(z, x)] : fz    
end



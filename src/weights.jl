# See Floater and Hormann 2007 Eqn. (18) (Actually, an unnumbered equation
# that is 3 equations later)
# UPDATE: used formula (18) from Hormann's "Barycentric Interpolation"
# part is that the interpolation seems to work OK for either. (We do know
# that the sign of the weights must alternate)
""" 
    floater_weights(x, d=0)  

Calculate the Floater-Hormann weights for grid x and mixing degree d.  For
many applications d = 3 or 4 works well.

    if we had a regular grid the values would be...
    x: points at which we have function values
    d: mixing degree        weights
    d = 0                 1,1,...,1,1  
    d = 1               1,2,2,...,2,2,1
    d = 2             1,3,4,4,...,4,4,3,1
    d = 3           1,4,7,8,8,...,8,8,7,4,1
    d = 4     1,5,11,15,16,16,...,16,16,15,11,5,1

But note, if you do have a regularly spaced grid, it is much better to call
the floater_weights(n, order)
"""
function floater_weights(x::Array{T, 1}, d=0) where {T}
    n = length(x)
    d > n && error("Mixing coefficient must be less than node set size")
    w = zeros(eltype(x), n)
    @inbounds for k in 1:n
        ws = 0.0
        for i in max(1, k-d):min(n-d, k)
            lwp = 0.0
            for j in i:min(n, i+d)
                j == k && continue
                lwp += log(abs(x[k] - x[j]))
            end
            ws += 1.0 / exp(lwp)
        end
        w[k] = k % 2 == 0 ? ws : -ws
    end
    return w ./ norm(w, Inf)
end

# Calculate the FH weights for a uniform grid of size n.  We don't
# worry about the actual grid size, since these weights are intended
# for use in the barycentric (2nd kind) interpolation formula.
# This is from the last equation (unnumbered) in Section 4 of the FH paper
function floater_weights(n::Int, d=0)
    w = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        s = 0.0
        for i in max(1, k-d):min(n-d, k)
            s += binomial(d, k-i)
        end
        w[k] = (k-d) % 2 == 0 ? s : -s
    end
    return w
end

# Recognize a step range and handle it efficiently
function floater_weights(x::T, d=0) where {T <:StepRangeLen}
    return floater_weights(length(x), d)
end


"""Return the n + 1 chebychev points on the [-1, 1] interval

  A shifted sin is used to make sure the points are symmetric
"""
function chebpts(n)
    m = n - 1
    x = sinpi.([-m:2:m;] ./ 2m)
    return x
end

function chebwts(n)
    w = ones(n)
    w[2:2:end] .= -1.0
    w[1] *= 0.5
    w[end] *= 0.5
    return w
end


# This is faster.  The exp(sum(log.(abs.(...))) version is ~30X slower do n = 10000
"""The weights for Barycentric Lagrange interpolation given the x grid of evaluations"""
function lagrange_weights(x::Vector{T}) where {T}
    n = length(x)
    w = ones(T, n)
    t = copy(w)
    @inbounds for i ∈ 1:n
        t .= x[i] .- x
        t[i] = one(eltype(x))
        w[i] = one(eltype(x)) / prod(t)
        # Here is the "safe" version
        # w[i] = one(eltype(x)) / exp(sum(log.(abs.(t))))
        # neg = count(t .< zero(eltype(x))) % 2 == 1
        # w[i] = neg ? -w[i] : w[i]
    end
    # Normalize weights
    return  w ./ norm(w, Inf)
end


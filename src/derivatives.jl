
# From Baltensperger and Trummer (eqns 3 - 5) we have the following.
function D1(r::T) where {T <: BRInterp}
    n = length(r.x)
    D = zeros(eltype(r.x), n, n)
    @inbounds for j ∈ 1:n
        for i ∈ 1:n
            if i != j
                D[i, j] = r.w[j] / (r.w[i] * (r.x[i] - r.x[j]))
            end
        end
        D[j, j] = - r.x[j] / (2*(one(eltype(r.x)) - r.x[j]^2))
    end
    D[1, 1] = (2*(n-1)^2 + 1) / 6
    D[end, end] = - D[1,1]
    return D
end


function differ(r::T) where {T <: FHInterp}
    n = length(r.x)

    # First the polynomial contribution
    D = D1(r)
    df = D * r.f

    # Now the rational "correction"
    dd2 = divided_difference(r.x ,r.f) ^ 2
    t = zeros(eltype(r.f), n)
    @inbounds for j ∈ 1:n
        t .= r.x[j] .- r.x
        idx = [1:j-1; j+1;n;]
        @views xv = r.x[idx]
        @views fv = r.f[idx]
        t = t ./ divided_difference(xv , fv)
        t[j] = one(eltype(r.x))
        df[j] += dd2 * prod(t)
    end
    
    dr = FHInterp(r.x, df, r.w, r.order)
end

function differ(r::T) where {T <: AAAapprox}
    D = D1(r)
    df = D * r.f
    dr = AAAapprox(r.x, df, r.w, zeros(eltype(r.x), 1))
end        

function divided_difference(x, y)
    n = length(x)
    s = zero(eltype(x))
    t = zeros(eltype(x), n)
    @inbounds for j ∈ 1:n
        t .= x[j] .- x
        t[j] = one(eltype(x))
        s += y[j] / prod(t)
    end
    return s
end

# From Proposition 12 of Schneider and Werner.
# Obviously, the default is for the first derivative.
# Also note that this works for both FHInterp and AAAapprox representations of
# rational polynomials
function deriv(a::T, z; m=1) where T <: BRInterp
    ET = eltype(a.x)
    V = a.w ./ (z .- a.x)
    γ = V ./ sum(V)
    δ = copy(a.f)
    ϕ = ET(0)
    for k in 0:m
        ϕ = γ' * δ
        δ .= (δ .- ϕ) ./ (a.x .- z)
    end
    r = factorial(m) * ϕ    
    isfinite(r) && return r

    j = findfirst(==(z), a.x)
    if j !== nothing
        γ = deleteat!(-a.w ./ a.w[j], j)
        δ = deleteat!((a.f .- a.f[j]) ./ (a.x .- a.x[j]), j)
        xx = copy(a.x)
        deleteat!(xx, j)
        for k in 1:m
            ϕ = γ' * δ            
            δ .= (δ .- ϕ) ./ (xx .- a.x[j])
        end
        return factorial(m) * ϕ    
    end
    error("Could not calculate derivative")
end

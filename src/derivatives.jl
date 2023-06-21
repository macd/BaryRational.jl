
# Evaluate the first derivative of interp object a at z.
function deriv(a::T, z) where T <: BRInterp
    V = a.w ./ (z .- a.x)
    denom = sum(V)
    U = (a(z) .- a.f) ./ (z .- a.x)
    numer = sum(V .* U)
    r = numer / denom
    isfinite(r) && return r

    i = findfirst(==(z), a.x)
    if i !== nothing
        s = 0.0
        @inbounds for j in eachindex(a.x)
            if j != i
                s += a.w[j] * ((a.f[j] - a.f[i]) / (a.x[j] - a.x[i]))
            end
        end
        return - s / a.w[i]
    end
    error("Could not calculate derivative")
end

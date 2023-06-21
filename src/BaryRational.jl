module BaryRational
export FHInterp, bary, chebpts, floater_weights, lagrange_weights
export aaa, prz, deriv

using LinearAlgebra
using GenericLinearAlgebra
using Statistics
using SparseArrays

abstract type BRInterp <: Function end

struct FHInterp{T <: AbstractFloat} <: BRInterp
    x::Vector{T}
    f::Vector{T}
    w::Vector{T}
    order::Int
end

include("weights.jl")
include("bary.jl")
include("aaa.jl")
include("derivatives.jl")

function FHInterp(x::AbstractVector{T}, f::Vector{T}; order::Int=0, grid=false) where {T <: AbstractFloat}
    if grid
        FHInterp(collect(x), f, floater_weights(length(x), T, d=order), order)
    else
        FHInterp(collect(x), f, floater_weights(x, d=order), order)
    end
end

(y::FHInterp)(z) = bary(z, y.f, y.x, y.w)


end

module BaryRational
export FHInterp, bary, chebpts, chebwts, floater_weights, lagrange_weights
export aaa, prz, deriv

using LinearAlgebra
using Statistics
using SparseArrays

abstract type BRInterp <: Function end

struct FHInterp{T <: AbstractArray} <: BRInterp
    x::T
    f::T
    w::T
    order::Int
end

include("weights.jl")
include("bary.jl")
include("aaa.jl")
include("derivatives.jl")

function FHInterp(x::Vector{T}, f::Vector{T}; order::Int=0, grid=false) where {T}
    if grid
        FHInterp(collect(x), f, floater_weights(length(x), order), order)
    else
        FHInterp(collect(x), f, floater_weights(x, order), order)
    end
end

(y::FHInterp)(z) = bary(z, y.f, y.x, y.w)


end

module BaryRational
export FHInterp, bary, chebpts, floater_weights, lagrange_weights
export aaa, prz, deriv

using GenericLinearAlgebra
using GenericSchur
using LinearAlgebra
using Printf
using SparseArrays
using Statistics

abstract type BRInterp <: Function end

struct FHInterp{T <: AbstractArray} <: BRInterp
    x::T
    f::T
    w::T
    order::Int
end

cplxord(t) = (real(t), imag(t))

include("weights.jl")
include("bary.jl")
include("aaa.jl")
include("derivatives.jl")

function FHInterp(x::AbstractVector{T}, f::AbstractVector{T}; order::Int=0, grid=false) where {T}
    if grid
        FHInterp(collect(x), f, floater_weights(length(x), T, d=order), order)
    else
        FHInterp(collect(x), f, floater_weights(x, d=order), order)
    end
end

(y::FHInterp)(z) = bary(z, y.f, y.x, y.w)


end

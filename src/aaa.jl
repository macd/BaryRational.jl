# AAA algorithm from the paper "The AAA Alorithm for Rational Approximation"
# by Y. Nakatsukasa, O. Sete, and L.N. Trefethen, SIAM Journal on Scientific
# Computing, 2018
using Printf

# so x can be real while f and w are complex
struct AAAapprox{T <: AbstractArray, W <: AbstractArray} <: BRInterp
    x::T
    f::W
    w::W
    errvec   # Interesting, but consider removing?
end

# NB: we don't sort errvec because that is just a history of the errors per
#     iteration of the algorithm
function Base.sort!(a::AAAapprox)
    ord = eltype(a.x) <: Complex ? cplxord : identity
    perm = sortperm(a.x, by=ord)
    permute!.((a.x, a.f, a.w), (perm,))
end

# In this version zz can be a scalar or a vector. BUT: do not make the mistake
# of broadcasting (ie a.(zz)) when zz is a vector. Although it gives correct
# results, it is much slower than just a(zz). We should really use bary rather
# than reval, but the tests are set up for this...
#(a::AAAapprox)(zz) = reval(zz, a.x, a.f, a.w)

# In this version zz is only ever a scalar, which means we must broadcast over
# the argument if it is a vector, BUT it looks to be much faster than reval.
# TODO: revisit
(a::AAAapprox)(zz) = bary(zz, a)
(a::AAAapprox)(zz::T) where {T <: AbstractVector} = bary.(zz, a)


function compute_weights(m, J, A::S) where {T, S <: AbstractMatrix{T}}
    if length(J) >= m                      # The usual tall-skinny case
        # Notice that A[J, :] selects only the non-support points and it is
        # those points that we will use for a least squares fit.
        G = svd(A[J, :])                   # Reduced SVD (the default)
        s = G.S
        mm = findall(==(minimum(s)), s)    # Treat case of multiple min sing val
        nm = length(mm)
        w = G.V[:, mm] * (ones(T, nm) ./ sqrt(T(nm))) # Aim for non-sparse wt vector
    elseif length(J) >= 1
        V = nullspace(A[J, :])             # Fewer rows than columns
        nm = size(V, 2)
        w = V * ones(T, nm) ./ sqrt(T(nm))   # Aim for non-sparse wt vector
    else
        w = ones(T, m) ./ sqrt(T(m))         # No rows at all
    end
    return w
end


"""
    aaa(Z, F; tol=1e-13, mmax=150, verbose=false, clean=1, do_sort=true) -> r::AAAapprox

Computes the rational approximation of data `F` on set `Z` using the AAA algorithm.

# Arguments
- `Z`: Vector of sample points.
- `F`: Vector of values at the points in `Z`.
- `tol`: Relative tolerance.
- `mmax`: Degree of numerator and denominator is at most `(mmax-1, mmax-1)`.
- `verbose`: If `true`, prints detailed information during computation.
- `clean`: If `true`, detects and removes Froissart doublets.
- `do_sort`: If `true` sorts the values of `Z` (and correspondingly `F`) in ascending order.

# Returns
- `r::AAAapprox`: A  struct representing the approximant that called as a function. The
   struct has fields, `z, f, w` = vectors of support points, function values, weights and
   `errvec` = vector of errors at each step.

Note 1. Changes from the MATLAB version include: Switched order of `Z` and `F` in the
function signature. Added `verbose` and `clean` boolean flags. Poles, residues, and zeros
(`pol`, `res`, `zer`) are calculated on demand by calling `prz(z::AAAapprox)`.

Note 2. The code (more or less) works with `BigFloat`. Since `prz` has not been made
generic, when using `BigFloat`, set `clean=false`. Specify `tol` as a tiny `BigFloat`
value explicitly, as default tolerances may not be sufficient.

# Example
```julia
    using BaryRational
    xrat = [-1//1:1//100:1//1;];
    xbig = BigFloat.(xrat);
    fbig = sin.(xbig);
    sf = aaa(xbig, fbig, verbose=true, clean=false, tol=BigFloat(1/10^40));

    julia @v1.10> sin(BigFloat(-1//3))
    -0.3271946967961522441733440852676206060643014068937597915900562770705763744817618

    julia @v1.10> sf(BigFloat(-1//3))
    -0.3271946967961522441733440852676206060643014068937597915900562770705763744817662
```
"""
function aaa(Z, F; tol=1e-13, mmax=150,
             verbose=false, clean=1, do_sort=true)
    U = eltype(Z)
    S = eltype(F)
    # filter out any NaN's or Inf's in the input
    keep = isfinite.(F)
    F = F[keep]
    Z = Z[keep]

    # Remove repeated elements of Z and corresponding elements of F
    ii = unique(i -> Z[i], eachindex(Z))
    Z = Z[ii]
    F = F[ii]

    M = length(Z)                    # number of sample points
    mmax = min(M, mmax)              # max number of support points

    abstol = tol * norm(F, Inf)
    verbose && println("\nabstol: ", abstol)

    F, Z = promote(F, Z)
    T = promote_type(S, U)

    J = [1:M;]
    z = T[]                          # support points
    f = T[]                          # function values at support points
    w = T[]                          # the weights
    A  = Matrix{T}(undef, M, 0)
    C  = Matrix{T}(undef, M, 0)

    errvec = T[]
    R = fill(mean(F), size(F))
    m = 1
    jtrunc = Int[]
    @inbounds for outer m in 1:mmax
        j = argmax(abs.(F .- R))     # select next support point
        push!(jtrunc, j)         # save index in case we need to truncate later
        push!(z, Z[j])
        push!(f, F[j])
        deleteat!(J, findfirst(isequal(j), J)) # update index vector

        # next column of Cauchy matrix
        C = hcat(C, T(1) ./ (Z .- Z[j]))

        # Loewner matrix
        A = hcat(A, (F .- f[end]) .* C[:, end])

        w = compute_weights(m, J, A)

        # Don't use the zero weights when calculating the approximation at the
        # support points
        i0 = findall(!=(T(0)), w)
        N = C[:, i0] * (w[i0] .* f[i0])       # numerator
        D = C[:, i0] *  w[i0]

        # Use the rational approximation at the remaining non support
        # points so we can measure the error.
        R .= F
        R[J] .= N[J] ./ D[J]

        err = norm(F - R, Inf)
        verbose && println("Iteration: ", m, "  err: ", err)
        push!(errvec, err)                    # max error at sample points
        err <= abstol && break                # stop if converged
    end

    # If we've gone to max iters, then it is possible that the best
    # approximation is at a smaller vector size. If so, truncate the
    # approximation which will give us a better approximation that is
    # faster to compute. Note that we must truncate A, reset J, and then
    # recompute the weights for this smaller size.
    if m == mmax
        idx = argmin(i -> real(errvec[i]), eachindex(errvec))
        if idx != mmax # if min error is at mmax, do nothing
            verbose && println("Hit max iters. Truncating approximation at $idx.")
            J = deleteat!([1:M;], sort(jtrunc))
            w = compute_weights(idx, J, @views(A[:, 1:idx]))
            for v in (z, f, errvec)
                deleteat!(v, idx+1:mmax)
            end
        end
    end

    # Remove the support points with zero weight.
    izero = findall(==(T(0)), w)
    for v in (z, f, w, errvec)
        deleteat!(v, izero)
    end

    r = AAAapprox(z, f, w, errvec)

    # We must sort if we plan on using bary rather than reval.
    do_sort && sort!(r)

    # Remove Froissart doublets if desired.  We do this in place.
    # TODO: use clean_tol instead of 1e-13
    if clean == 1
        cleanup!(r, Z, F; verbose=verbose)
    elseif clean == 2
        r = cleanup2!(r, Z, F; verbose=verbose)
    elseif clean == 3
        old_cleanup!(r, Z, F; verbose=verbose)
    end

    return r
end

function prz(r::AAAapprox)
    return prz(r.x, r.f, r.w)
end


# We now use the Schur decompostion so we can use GenericShur which means that
# we can now use prz with BigFloat,
function prz(z, f, w)
    T = eltype(z)
    m = length(w)
    B = diagm(ones(T, m+1))
    B[1, 1] = T(0)
    E = [T(0)  transpose(w); ones(T, m) diagm(z)]
    sp = schur(Complex.(E), Complex.(B))
    pol = sp.values[isfinite.(sp.values)]

    # Compute residues as quotient of analytic functions
    N =   (T(1) ./ (pol .- transpose(z))) * (f .* w)
    D = -((T(1) ./ (pol .- transpose(z))) .^ 2) * w
    res = N ./ D

    E = [T(0) transpose(w .* f); ones(T, m) diagm(z)];
    sz = schur(Complex.(E), Complex.(B))
    zer = sz.values[isfinite.(sz.values)]

    pol, res, zer
end

# This is just the barycentric interpolation formula in matrix form.
# It is doing the same calculation as bary(...) only it does not need
# to be broadcasted over zz. (and in fact, you should not do so)
function reval(zz, z, f, w)
    # evaluate r at zz
    zv = size(zz) == () ? [zz] : vec(zz)
    CC = 1.0 ./ (zv .- transpose(z))         # Cauchy matrix
    r = (CC * (w .* f)) ./ (CC * w)          # AAA approx as vector
    r[isinf.(zv)] .= sum(f .* w) ./ sum(w)

    ii = findall(isnan.(r))               # find values NaN = Inf/Inf if any
    @inbounds for j in ii
        # Wow, linear search, but only if a NaN happens
        v = findfirst(==(zv[j]), z)
        if !isnan(zv[j]) && (v !== nothing)
            r[j] = f[v]  # force interpolation there
        end
    end
    r = size(zz) == () ? r[1] : reshape(r, size(zz))  # the AAA approximation
end


# This modifies the rational approximant r, if necessary. Updated July 2023 to
# more or less match the Chebfun version.
function cleanup!(r, Zp::AbstractVector{T}, Fp::AbstractVector{T};
                  verbose=false, cleanup_tol=1e-13) where {T}
    z, f, w = copy(r.x), copy(r.f), copy(r.w)
    pol, res, zer = prz(z, f, w)

    # Don't modify the original input vectors
    Z = copy(Zp)
    F = copy(Fp)

    # find negligible residues
    if any(F .!= T(0))
        geometric_mean_of_absF = exp(mean(log.(abs.(F[F .!= T(0)]))))
    else
        geometric_mean_of_absF = T(0)
    end

    P = length(pol)
    zdistances = Vector{typeof(abs(T(0)))}(undef, P)
    for j in 1:P
        zdistances[j] = norm(pol[j] .- Z, -Inf)
    end
    ii = findall(abs.(res) ./ zdistances .< cleanup_tol * geometric_mean_of_absF)

    ni = length(ii)
    ni == 0 && return
    sn = ni == 1 ? "" : "s"
    verbose && println("$ni Froissart doublet$sn. Number of residues = ", length(res))

    # For each spurious pole find and remove closest support point:
    @inbounds for j = 1:ni
        azp = abs.(z .- pol[ii[j]] )
        _, jj = findmin(azp)
        deleteat!(z, jj)    # remove nearest support points
        deleteat!(f, jj)
    end

    # Remove support points z from sample set:
    @inbounds for zs in z
        idx = findfirst(==(zs), Z)
        deleteat!(F, idx)
        deleteat!(Z, idx)
    end
    m = length(z)
    M = length(Z)
    SF = spdiagm(M, M, 0 => F)
    Sf = diagm(f)
    C = 1 ./ (Z .- transpose(z))
    A = SF*C - C*Sf
    G = svd(A)
    w = G.V[:, m]

    if length(z) < length(r.x)
        resize!.((r.x, r.f, r.w), length(z))
        r.x .= z
        r.f .= f
        r.w .= w
    end
    return r
end

# The old ways are sometimes the good ways... This was coded from the original
# AAA paper. Only calculate the updated z, f, and w
function old_cleanup!(r, Zp::AbstractVector{T}, Fp::AbstractVector{T};
                  verbose=false, cleanup_tol=1e-13) where {T}
    z, f, w = r.x, r.f, r.w
    pol, res, zer = prz(z, f, w)
    Z = copy(Zp)
    F = copy(Fp)
    m = length(z)
    M = length(Z)
    ii = findall(abs.(res) .< 1e-13)  # find negligible residues
    ni = length(ii)
    ni == 0 && return
    println("$ni Froissart doublets. Number of residues = ", length(res))

    # For each spurious pole find and remove closest support point:
    @inbounds for j = 1:ni
        azp = abs.(z .- pol[ii[j]] )
        jj = findall(isequal(minimum(azp)), azp)
        deleteat!(z, jj)    # remove nearest support points
        deleteat!(f, jj)
    end

    # Remove support points z from sample set:
    @inbounds for j = 1:length(z)
        jj = findall(isequal(z[j]), Z)
        deleteat!(F, jj)
        deleteat!(Z, jj)
    end
    m = m - length(ii)
    SF = spdiagm(M-m, M-m, 0 => F)
    Sf = diagm(f)
    C = 1 ./ (Z .- transpose(z))
    A = SF*C - C*Sf
    G = svd(A)
    ww = G.V[:, m]
    deleteat!(r.w, 1:(length(r.w) - length(ww)))
    r.w .= ww
    return nothing
end


# Alternative cleanup procedure to remove spurious pole-zero pairs.
# This considers pole-zero distances. Ported from Chebfun July 2023
# Note that we pass in the original sample set instead of sticking it
# on the struct.
function cleanup2!(r, Zp::AbstractVector{T}, Fp::AbstractVector{T};
                   cleanup_tol=1//10^13, verbose=false) where {T}
    z = copy(r.x); f = copy(r.f); w = copy(r.w)
    pol, res, zer = prz(z, f, w)
    FT = typeof(abs(T(0)))
    cleanup_tol = FT(cleanup_tol)

    niter = 0
    while true
        niter = niter + 1
        Z = copy(Zp)
        F = copy(Fp)
        ii = Int[]
        @inbounds for jj in 1:length(pol)
            dz = isempty(zer) ? FT(10^100) : minimum(abs.(zer .- pol[jj]))
            dS = abs.(Z .- pol[jj])
            ds = minimum(dS)
            if any(F .!= T(0))
                q = 4 * pi * abs.(F) .* dS
                Q = mean(q)                # Arithmetic mean
            else
                Q = T(0)
            end
            R = 8 * cleanup_tol * Q / 4pi   # Equivalent residue value

            # Conditions to expunge poles
            # Expunge if either minimum distance is zero
            if ds == FT(0) || dz == FT(0)
                push!(ii, jj)
            # Expunge if Z is a real interval
            elseif isreal(Z) && abs(imag(pol[jj])) < eps(FT) &&
                (real(pol[jj]) >= minimum(Z)) && (real(pol[jj]) <= maximum(Z))
                push!(ii, jj)
            # Expunge if Z is the unit disk
            elseif all(abs.(Z) == T(1)) && abs(abs(pol[jj]) - 1) < eps(FT)
                push!(ii, jj)
            # Expunge if distance to closest zero is undetectable
            elseif (dz/ds) < FT(1) && dz < max(cleanup_tol^2, eps(FT))
                push!(ii, jj)
            # Expunge if a nearby zero exists and residue is below the
            # equivalent value R. Two choices for real and complex F
            elseif (dz/ds) < sqrt(cleanup_tol)
                if  T <: AbstractFloat && abs(real(res[jj])) < R
                    push!(ii, jj)
                elseif abs(res[jj]) < R
                    push!(ii, jj)
                end
            end
        end
        unique!(ii)

        ni = length(ii)
        if ni == 0
            # Nothing to do.
            break
        else
            ns = ni == 1 ? "" : "s"
            verbose && println("AAA:Froissart:  ",
                               "$ni Froissart doublet$ns, niter = ", niter)
        end

        # For each spurious pole find and remove closest support point:
        for j = 1:ni
            azp = abs.(z .- pol[ii[j]])
            _, jj = findmin(azp)
            verbose && println("Removing point at ", z[jj])
            # Remove support point(s):
            deleteat!(z, jj)
            deleteat!(f, jj)
        end

        # Remove support points z from sample set:
        @inbounds for zs in z
            idx = findfirst(==(zs), Z)
            deleteat!(F, idx)
            deleteat!(Z, idx)
        end
        m = length(z)
        M = length(Z)

        # Build Loewner matrix:
        SF = spdiagm(M, M, 0 => F)
        Sf = diagm(f)
        C = 1 ./ (Z .- transpose(z))   # Cauchy matrix.
        A = SF * C - C * Sf            # Loewner matrix.

        # Solve least-squares problem to obtain weights:
        G = svd(A)
        @views w = G.V[:, m]

        # Compute poles, residues and zeros for next round.
        pol, res, zer = prz(z, f, w)
    end   # End of while loop

    # reconstruct r
    if length(z) < length(r.x)
        resize!.((r.x, r.f, r.w), length(z))
        r.x .= z
        r.f .= f
        r.w .= w
    end
    return r
end

# Borrowed from https://github.com/complexvariables/RationalFunctionApproximation.jl

#####
##### Adaptive AAA on [-1, 1] only
#####

# refinement in parameter space
function refine(t, N)
    x = sort(t)
    Δx = diff(x)
    d = eltype(x).((1:N) / (N+1))
    return vec( x[1:end-1] .+ (d' .* Δx) )
end

function aaa(
    f::Function, ::Type{T}=Float64;
    mmax=150, tol=1000*eps(T),
    refinement=3, lookahead=10, stats=false
    ) where {T <: AbstractFloat}
    CT = Complex{T}
    # arrays for tracking convergence progress
    err, nbad = T[], Int[]
    nodes, vals, pol, weights = Vector{T}[], Vector{CT}[], Vector{CT}[], Vector{CT}[]

    S = [-one(T), one(T)]                       # initial nodes
    fS = f.(S)
    besterr, bestm = Inf, NaN
    while true                                  # main loop
        m = length(S)
        push!(nodes, copy(S))
        X = refine(S, max(refinement, ceil(16-m)))    # test points
        fX = f.(X)
        push!(vals, copy(fS))
        C = [ 1/(x-s) for x in X, s in S ]
        L = [a-b for a in fX, b in fS] .* C
        _, _, V = svd(L)
        w = V[:,end]
        push!(weights, w)
        R = (C*(w.*fS)) ./ (C*w)                # values of the rational interpolant
        push!(err, norm(fX - R, Inf) )

        #zp =  poles(Barycentric(S, fS, w))
        zp, _, _ = prz(S, fS, w)
        push!(pol, zp)
        I = (imag(zp).==0) .& (abs.(zp).<=1)    # bad poles indicator
        push!(nbad, sum(I))
        # If valid and the best yet, save it:
        if (last(nbad) == 0) && (last(err) < besterr)
            besterr, bestm = last(err), m
        end

        fmax = max( norm(fS, Inf), norm(fX, Inf) )     # scale of f
        # Check stopping:
        if (besterr <= tol*fmax) ||                             # goal met
            (m == mmax + 1) ||                                  # max degree reached
            ((m - bestm >= lookahead) && (besterr < 1e-2*fmax)) # stagnation
            break
        end

        # We're continuing the iteration, so add the worst test point to the nodes:
        _, j = findmax(abs, fX - R)
        push!(S, X[j])
        push!(fS, fX[j])
    end

    # Use the best result found:
    S, y, w = nodes[bestm-1], vals[bestm-1], weights[bestm-1]
    idx = sortperm(S)
    x, y, w = S[idx], y[idx], w[idx]
    if isreal(w) && isreal(y)
        y, w = real(y), real(w)
    end

    # TODO: Consider borrowing the stats stuff too.
    # if stats
    #     if isreal(w) && isreal(y)
    #         weights = real.(weights)
    #         vals = real.(vals)
    #     end
    #     st = ConvergenceStats(bestm-1, err, nbad, nodes, vals, weights, pol)
    #     r = AAAapprox(x, y, w, [])
    # else
    #     r = AAAapprox(x, y, w)
    # end

    r = AAAapprox(x, y, w, [])
    return r
end

# AAA algorithm from the paper "The AAA Alorithm for Rational Approximation"
# by Y. Nakatsukasa, O. Sete, and L.N. Trefethen, SIAM Journal on Scientific Computing
# 2018

struct AAAapprox{T <: AbstractArray} <: BRInterp
    x::T
    f::T
    w::T
    errvec::T
end

# In this version zz can be a scalar or a vector
(a::AAAapprox)(zz) = reval(zz, a.x, a.f, a.w)

# In this version zz is only ever a scalar
#(a::AAAapprox)(zz) = bary(zz, a)

# Handle function inputs as well.  Seems unnecessary, consider removing.
function aaa(Z::AbstractArray{T,1}, F::S;  tol=1e-13, mmax=100, verbose=false, clean=true, do_sort=false) where {T, S<:Function}
    aaa(Z, F.(Z), tol=tol, mmax=mmax, verbose=verbose, clean=clean, do_sort=do_sort)
end

"""
    aaa(Z, F; tol=1e-13, mmax=100, verbose=false, clean=true, do_sort=false)

Computes a rational approximant of the data F on the set Z using the AAA algorithm.

# Arguments 
- `Z`: The vector of sample points.
- `F`: The vector of data values corresponding to the points in `Z`. `F` can also be a function, in which case `F.(Z)` is used.

# Keyword Arguments 
- `tol=1e-13`: The relative tolerance.
- `mmax=100`: Sets the maximum type of the rational approximant to `(mmax-1, mmax-1)`.
- `verbose=false`: Prints info while calculating if `true`.
- `clean=true`: Detects and remove Froissart doublets if `true`.
- `do_sort=false`: Sorts the support points if `true`. Does not work if the sample data are complex. 

# Output 
The returned value is an approximant `r::AAAapprox` that could be called e.g. as `r(z)`, returning the approximant at `z`.
See the docs for examples.

!!! note "Changes from MATLAB version"

    - The order of `Z` and `F` are changed in the function signature.
    - The `verbose` and `clean` boolean flags are added.
    - The vectors of poles, residues, and zeros are now only calculated on demand by calling `prz(r::AAAapprox)`.

!!! note "BigFloats"

    This does (more or less) work with `BigFloat`s. 
    Caveats: since `prz` has not yet been made generic, you must set `clean=false`. 
    Also, you set `tol` to a tiny `BigFloat` value rather than use the defaults. For example:

    ```julia-repl 
    julia> using BaryRational
    julia> xrat = -(1//1):(1//100):(1//1);
    julia> xbig = BigFloat.(xrat);
    julia> fbig = sin.(xbig);
    julia> sf = aaa(xbig, fbig, clean = false, tol = BigFloat(1//10^40));
    julia> sin(BigFloat(-1//3))
    -0.3271946967961522441733440852676206060643014068937597915900562770705763744817618
    
    julia> sf(BigFloat(-1//3))
    -0.3271946967961522441733440852676206060643014068937597915900562770705763744817662
    ```

"""
function aaa(Z::AbstractArray{U,1}, F::AbstractArray{S,1}; tol=1e-13, mmax=100,
             verbose=false, clean=true, do_sort=false) where {S, U}
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
    @inbounds for outer m in 1:mmax
        j = argmax(abs.(F .- R))               # select next support point
        push!(z, Z[j])
        push!(f, F[j])
        deleteat!(J, findfirst(isequal(j), J)) # update index vector

        # next column of Cauchy matrix
        C = hcat(C, T(1) ./ (Z .- Z[j]))

        # Loewner matrix
        A = hcat(A, (F .- f[end]) .* C[:, end])

        # Compute weights:
        if length(J) >= m                      # The usual tall-skinny case
            G = svd(A[J, :])                   # Reduced SVD (the default)
            s = G.S
            mm = findall(==(minimum(s)), s)    # Treat case of multiple min sing val
            nm = length(mm)
            w = G.V[:, mm] * (ones(T, nm) ./ sqrt(nm)) # Aim for non-sparse wt vector
        elseif length(J) >= 1
            V = nullspace(A[J, :])             # Fewer rows than columns
            nm = size(V, 2)                    
            w = V * ones(T, nm) ./ sqrt(nm)   # Aim for non-sparse wt vector
        else
            w = ones(T, m) ./ sqrt(m)         # No rows at all (needed for Octave)
        end

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
        push!(errvec, err)                    # max error at sample point
        err <= abstol && break                # stop if converged
    end

    # If we've gone to max iters, then it is possible that the best 
    # approximation is at a smaller vector size. If so, truncate the
    # approximation which will give us a better approximation that is
    # faster to compute.
    if m == mmax
        verbose && println("Hit max iters. Truncating approximation.")
        idx = argmin(errvec)
        deleteat!(z, idx+1:mmax)
        deleteat!(f, idx+1:mmax)
        deleteat!(w, idx+1:mmax)
        deleteat!(errvec, idx+1:mmax)
    end

    # Remove the support points with zero weight.
    izero = findall(==(T(0)), w)
    deleteat!(z, izero)
    deleteat!(f, izero)
    deleteat!(w, izero)
    deleteat!(errvec, izero)
    
    # We must sort if we plan on using bary rather than reval, _but_ this
    # will not work when z is complex
    if do_sort
        perm = sortperm(z)
        z .= z[perm]
        f .= f[perm]
        w .= w[perm]
        errvec .= errvec[perm]
    end
    r = AAAapprox(z, f, w, errvec)

    # Remove Froissart doublets if desired.  We do this in place, but must
    # skip this step (for now, set clean=false) if we are using BigFloats
    if clean
        pol, res, zer = prz(r)            # poles, residues, and zeros
        ii = findall(abs.(res) .< 1e-13)  # find negligible residues
        length(ii) != 0 && cleanup!(r, pol, res, zer, Z, F)
    end
    
    return r
end


# NB: So this is not (yet) set up to be generic. If trying to use aaa with
# BigFloat's be sure to set clean=false.
"""
    prz(r::AAAapprox)

Return the poles, residues, and zeros of the AAA approximation r. 

!!! warning 

    If you are trying to use `aaa` with `BigFloat`s, make sure you have called 
    `aaa` with `clean=false` as, currently, `prz` is not generic.
"""
function prz(r::AAAapprox)
    z, f, w = r.x, r.f, r.w
    T = eltype(z)
    m = length(w)
    B = diagm(ones(T, m+1))
    B[1, 1] = T(0)
    E = [T(0)  transpose(w); ones(T, m) diagm(z)]
    pol, _ = eigen(E, B)
    pol = pol[isfinite.(pol)] 
    dz = T(1//100000) * exp.(2im*pi*[1:4;]/4)
    
    # residues
    res = r(pol .+ transpose(dz)) * dz ./ 4 
        
    E = [0 transpose(w .* f); ones(T, m) diagm(z)]
    zer, _ = eigen(E, B)
    zer = zer[isfinite.(zer)]
    pol, res, zer
end

# This is just the barycentric interpolation formula in matrix form.
# It is doing the same calculation as bary(...) only it does not need
# to be broadcasted over zz. For a random vector xx of length 1000 between
# -1 and 1, with a function with 17 support points, I see the following:
#
# julia @v1.10> @btime dya = g(xx);
#   20.143 μs (15 allocations: 165.98 KiB)
#
# julia @v1.10> @btime dya = bary.(xx, g);
#   17.230 μs (6 allocations: 8.12 KiB)
#
# Also, this version allocates like crazy. Maybe set the default to bary?
#
function reval(zz, z, f, w)
    # evaluate r at zz
    zv = size(zz) == () ? [zz] : vec(zz)  
    CC = 1.0 ./ (zv .- transpose(z))         # Cauchy matrix
    r = (CC * (w .* f)) ./ (CC * w)          # AAA approx as vector
    r[isinf.(zv)] .= sum(f .* w) ./ sum(w)
    
    ii = findall(isnan.(r))               # find values NaN = Inf/Inf if any
    @inbounds for j in ii
        # Wow, linear search , but only if a NaN happens
        v = findfirst(==(zv[j]), z)
        if !isnan(zv[j]) && (v !== nothing)
            r[j] = f[v]  # force interpolation there
        end
    end
    r = size(zz) == () ? r[1] : reshape(r, size(zz))  # the AAA approximation
end


# Only calculate the updated z, f, and w
# FIXME: Change the hardcoded tolerance in this function
function cleanup!(r, pol, res, zer, Z, F)
    z, f, w = r.x, r.f, r.w
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


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
function aaa(Z::AbstractVector{T}, F::S;  tol=1e-13, mmax=100, verbose=false,
             clean=true) where {T, S<:Function}
    aaa(Z, F.(Z), tol=tol, mmax=mmax, verbose=verbose, clean=clean)
end

"""aaa  rational approximation of data F on set Z
        r = aaa(Z, F; tol, mmax, verbose, clean)

 Input: Z = vector of sample points
        F = vector of data values at the points in Z
        tol = relative tolerance tol, set to 1e-13 if omitted
        mmax: max type is (mmax-1, mmax-1), set to 100 if omitted
        verbose: print info while calculating default = false
        clean: detect and remove Froissart doublets default = true

 Output: r = an AAA approximant as a callable struct with fields
         z, f, w = vectors of support pts, function values, weights
         errvec = vector of errors at each step

 Note 1: Changes from matlab version:
         switched order of Z and F in function signature
         added verbose and clean boolean flags
         pol, res, zer = vectors of poles, residues, zeros are now only 
         calculated on demand by calling prz(z::AAAapprox)

 Note 2: This does (more or less) work with BigFloats. Caveats: since prz
         has not been made generic, you must set clean=false. Also, must
         set tol to a tiny BigFloat value rather than use the defaults.

    using BaryRational
    xrat = [-1//1:1//100:1//1;];
    xbig = BigFloat.(xrat);
    fbig = sin.(xbig);
    sf = aaa(xbig, fbig, verbose=true, clean=false, tol=BigFloat(1/10^40));

    julia @v1.10> sin(BigFloat(-1//3))
    -0.3271946967961522441733440852676206060643014068937597915900562770705763744817618

    julia @v1.10> sf(BigFloat(-1//3))
    -0.3271946967961522441733440852676206060643014068937597915900562770705763744817662
"""
function aaa(Z::AbstractVector{U}, F::AbstractVector{S}; tol=1e-13, mmax=100,
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
            # Notice that A[J, :] selects only the non-support points
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
            w = ones(T, m) ./ sqrt(m)         # No rows at all
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
        for v in (z, f, w, errvec)
            deleteat!(v, idx+1:mmax)
        end
    end

    # Remove the support points with zero weight.
    izero = findall(==(T(0)), w)
    for v in (z, f, w, errvec)
        deleteat!(v, izero)
    end
        
    # We must sort if we plan on using bary rather than reval, _but_ this
    # will not work when z is complex
    if do_sort
        perm = sortperm(z)
        for v in (z, f, w, errvec)
            permute!(v, perm)
        end
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
        # Wow, linear search, but only if a NaN happens
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

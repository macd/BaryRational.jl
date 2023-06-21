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

# Handle function inputs as well
function aaa(Z::AbstractArray{T,1}, F::S;  tol=1e-13, mmax=100, verbose=false, clean=true) where {T, S<:Function}
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
function aaa(Z::AbstractArray{T,1}, F::AbstractArray{S,1}; tol=1e-13, mmax=100,
             verbose=false, clean=true, do_sort=false) where {S, T}
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
    
    reltol = tol * norm(F, Inf)
    verbose && println("\nreltol: ", reltol)

    SF = spdiagm(M, M, 0 => F)       # left scaling matrix
    
    F, Z = promote(F, Z)
    P = promote_type(S, T)
    
    J = [1:M;]
    z = P[]                          # support points
    f = P[]                          # function values at support points
    C = P[]
    w = P[]
    
    errvec = P[]
    R = fill(mean(F), size(F))
    @inbounds for m = 1:mmax
        j = argmax(abs.(F .- R))               # select next support point
        push!(z, Z[j])
        push!(f, F[j])
        deleteat!(J, findfirst(isequal(j), J))   # update index vector

        # next column of Cauchy matrix
        C = isempty(C) ? reshape((1 ./ (Z .- Z[j])), (M,1)) : [C (1 ./ (Z .- Z[j]))]

        Sf = diagm(f)                         # right scaling matrix
        A = SF * C - C * Sf                   # Loewner matrix

        # There are times when we need the full decomposition here. We could
        # just always set full=true, but that slows down the test suite by
        # almost 3X
        mv, nv = size(A[J, :])
        G = svd(A[J, :], full=(mv < nv))
        
        w = G.V[:, m]
        
        N = C * (w .* f)                      # numerator 
        D = C * w
        R .= F
        R[J] .= N[J] ./ D[J]                  # rational approximation
        
        err = norm(F - R, Inf)
        verbose && println("Iteration: ", m, "  err: ", err)
        errvec = [errvec; err]                # max error at sample points
        err <= reltol && break                # stop if converged
    end

    # we must sort if we plan on using bary rather than reval
    if do_sort
        perm = sortperm(z)
        z .= z[perm]
        f .= f[perm]
        w .= w[perm]
    end
    r = AAAapprox(z, f, w, errvec)

    # remove Frois. doublets if desired.  We do this in place
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
    m = length(w)
    B = diagm(ones(m+1))
    B[1, 1] = 0.0
    E = [0.0  transpose(w); ones(m) diagm(z)]
    pol, _ = eigen(E, B)
    pol = pol[isfinite.(pol)] 
    dz = 1e-5 * exp.(2im*pi*[1:4;]/4)
    
    # residues
    res = r(pol .+ transpose(dz)) * dz ./ 4 
        
    E = [0 transpose(w .* f); ones(m) diagm(z)]
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
        if !isnan(zv[j]) && ((v = findfirst(==(zv[j]), z)) !== nothing)
            r[j] = f[v]  # force interpolation there
        end
    end
    r = size(zz) == () ? r[1] : reshape(r, size(zz))  # the AAA approximation
end


# Only calculate the updated z, f, and w
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


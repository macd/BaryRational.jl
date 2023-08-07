# This is a Julia translation of the code that appears in "AAA-least
# squares rational approximation and solution of Laplace problems" by Costa
# and Trefethen. If you just cut and paste their code into Matlab, it works
# as advertised. If you run it in Octave, it runs correctly (I assume) but it
# gives worse results. The results from the code here are very close to the 
# Matlab version. For instance, the max errors on the boundary only start 
# to differ in the 5th place. The default size of the plot here is too small,
# so to see all the details, expand it.

using LinearAlgebra
using BaryRational
using Base.MathConstants
using PyPlot

# Well, inpolygon from GeometricalPedicates only works for Convex polygons.
# Bummer. So here we just do successive testing on whether the point is in
# any of 2 regions that can be used to make the L shaped region.
# Just to be explicit, this works only for the L shaped region of the paper.
function inpolygonc(z::Complex{T}) where {T}
    x, y = reim(z)
    # Bottom rectangle
    (T(0) <= x <= T(2)) && (T(0) <= y <= T(1)) && return true
    # Top square
    (T(0) <= x <= T(1)) && (T(1) <= y <= T(2)) && return true
    return false
end

# "In MATLAB we use contructions like logspace(-14, 0, 300)' for a singularity
#  at one endpoint of [0, 1] and tanh(linspace(-16,16,600)') for singularities
#  at both endpoints of [−1, 1]."
#
# NOTE: I set this up to be generic so we can test different precisions.
#       The results of that testing resulted in a SingularException in the
#       calculation of A \ H. To get around that, for Float128, I just did
#       pinv(A) * B. For a tol=1e-12, it results in errors that all seem to be
#       an order of magnitude larger. Note also that in the paper, they talk 
#       about "the regularizing effects of least-squares solvers as realized 
#       in the MATLAB backslash command." There is more going on here, but 
#       the rabbit hole is already too deep. Time to declare victory and
#       move on.
function main(::Type{T}=Float64; tol=T(1//10^8)) where {T}
    # Set up
    plt.rcParams["text.usetex"] = true
    s = complex(tanh.(T.(range(-12//1, 12//1, length=300))))
    Z = Complex{T}.([ 1//1           .+          s;
                     (2//1 + 1im//2) .+ 1im//2 * s;
                     (3//2 + 1im//1) .+ 1//2   * s;
                     (1//1 + 3im//2) .+ 1im//2 * s;
                     (1//2 + 2im//1) .+ 1//2   * s;
                      1im//1         .+ 1im    * s])

    # These are the corners of the L shaped region
    w = Complex{T}.([0, 2, 2 + 1im, 1 + 1im, 1 + 2im, 2im])

    # We don't really need an anonymous function here since we don't reuse it.
    h = z -> real.(z) .^ 2
    H = h(Z)

    # Local AAA fits
    #tol = T(1//10^8)   # Interesting... when tol=1e12 we get the loopy poles
                        # and a lower accuracy
    pol_in  = Complex{T}[]
    pol_out = Complex{T}[]
    ax1 = subplot(221)
    for k = 1:6
        println("starting iteration: ", k)
        # this finds the points in Z that are closest to w[k] vs. than say,
        # any other w[j].
        ii = findall(abs.(Z .- w[k]) .== minimum(abs.(Z .- transpose(w)), dims=2))
        
        # This plots the sample points, colored by closest corner, so we 
        # can see which sample points are used by the local AAA 
        # approximation to determine the local poles.
        ax1.plot(real(Z[ii]), imag(Z[ii]), "*", markersize=3)
        
        g = aaa(Z[ii], H[ii], tol=tol, clean=0)
        polk, _, _ = prz(g)
        polk_in = polk[inpolygonc.(polk)]
        append!(pol_in, polk_in)
        polk_out = polk[.!inpolygonc.(polk)]
        append!(pol_out, polk_out)
    end
    ax1.axis("square")

    # Somehow, we calculate a different (from Costa & Trefethen) set of poles 
    # around the third corner. Here we see only 67 poles that bisect the
    # corner. In the Matlab code, with 98 poles, the poles have that funny
    # loopy shape. All the other corners appear to be roughly equivalent (at
    # least they result in the same number of poles)
    ax1.plot(real(pol_out), imag(pol_out), ".r", markersize=6)
    ax1.plot(real(pol_in),  imag(pol_in),  ".b", markersize=6)
    
    title("Local AAA poles")
    xlim(-.8, 2.8)
    ylim(-.8, 2.8)
    
    pol = transpose(pol_out)
    d = minimum(abs.(w .- pol), dims=1)
    println("pol: ", length(pol), "  d: ", length(d))

    # Solution
    Hes, P = VAorthog(Z, 20)
    Q = d ./ (Z .- pol)
    A = [real(P) real(Q) -imag(P) -imag(Q)]
    c = reshape(A \ H, :, 2) * Complex{T}.([1, 1im])

    # The following is needed for Float128 (and presumably BigFloat). See 
    # the comments above
    #c = reshape(pinv(A) * H, :, 2) * Complex{T}.([1, 1im])
    
    println("P: ", size(P), "  Q:", size(Q), "  A: ", size(A), "  c: ", size(c))
    F = [P Q] * c
    U = real(F)
    println("max error on boundary: ", norm(H - U, Inf))
    f = z -> reshape([VAeval(z[:], Hes)[1] d./(z[:] .- pol)]*c, size(z))
    u = z -> real(f(z))

    # Contour and error plots
    ax2 = subplot(222)
    ax2.axis("square")
    x, y = reim([w[:]; w[1]])
    ax2.plot(x, y)
    ax2.plot(real(pol), imag(pol), ".r", markersize=6)
    title("Laplace solution")
    xlim(-.8, 2.8)
    ylim(-.8, 2.8)

    N = 150
    x = LinRange(0.0, 2.0, N)
    xx = x' .* ones(N)
    yy = ones(N)' .* x
    zz = complex.(xx, yy)
    uu = u(zz)
    uu[.!inpolygonc.(zz)] .= NaN
    contour(x, x, uu, 20)

    ax3 = subplot(223)
    ax3.axis("square")
    ax3.semilogy(angle.(Z.-(.5+.5im)), abs.(U - H), ".k", markersize=3)
    xlim(-pi, pi)
    ylim(1e-12, 1e-4)
    grid(true, which="both")
    ax3.set_xticks([pi*(-1:.5:1);])
    ax3.set_xticklabels([raw"-$\pi$", raw"-$\pi$/2", "0", raw"$\pi$/2", raw"$\pi$"])
    title("Error against angle")
    nothing
end

# Vandermonde + Arnoldi orthogonalization
#   Input:
#      Z = column vector of sample points
#      n = degree of polynomial (>= 0)
#      Pol = array of vectors of poles (optional)
#   Output:
#      Hes = array of Hessenberg matrices (length 1+length(Pol))
#      R = matrix of basis vectors
function VAorthog(Z::AbstractVector{T}, n, Pol=T[]) where {T}
    M = length(Z)
    Q = ones(T, M, n+1)
    H = zeros(T, n+1, n)
    for k in 1:n
        q = Z .* Q[:, k]
        for j in 1:k
            H[j, k] = Q[:, j]' * q / M
            q .= q - H[j, k] * Q[:, j]
        end
        H[k+1, k] = norm(q) / sqrt(M)
        Q[:, k+1] .= q / H[k+1, k]
    end
    Hes = []  # yes Any[]. It corresponds to a "cell" array in Matlab
    push!(Hes, H)
    R = Q
    # Next orthogonalize the pole parts, if any
    TPol = copy(Pol)
    while !isempty(TPol)
        pol = popfirst!(TPol)
        np = length(pol)
        H = zeros(T, np, np-1)
        Q = ones(T, M, np+1)
        for k in 1:np
            q = Q[:, k] ./ (Z .- pol[k])
            for j in 1:k
                H[j, k] = Q[:, j]' * q ./ M
                q .= q - H[j, k] .* Q[:, j]
            end
            H[k+1, k] = norm(q) / sqrt(M)
            Q[:, k+1] .= q ./ H[k+1, k]
        end
        push!(Hes, H)
        R = hcat(R, Q[:, 2:end])
    end
    return Hes, R
end

# Vandermonde + Arnoldi basis construction
#   Input:
#      Z = column vector of sample points
#      Hes = array of Hessenberg matrices
#      Pol = array of vectors of poles, if any
#   Output:
#      R0 = matrix of basis vectors for functions
#      R1 = matrix of basis vectors
function VAeval(Z::AbstractVector{T}, Hes, Pol = T[]) where {T}
    M = length(Z)
    HE = deepcopy(Hes)
    H = popfirst!(HE)
    n = size(H, 2)
    Q = ones(T, M, n+1)
    D = zeros(T, M, n+1)

    for k = 1:n
        hkk = H[k+1, k]
        Q[:, k+1] .= (Z .* Q[:, k] .- Q[:, 1:k] * H[1:k, k]) / hkk
        D[:, k+1] .= (Z .* D[:, k] .- D[:, 1:k] * H[1:k, k] .+ Q[:, k]) / hkk
    end

    R0 = Q
    R1 = D
    # Next construct the pole parts of the basis, if any
    TPol = copy(Pol)
    while !isempty(TPol)
        pol = popfirst!(TPol)
        H = popfirst!(HE)
        np = length(pol)
        Q = ones(T, M, np)
        D = zeros(T, M, np)
        for k in 1:np
            Zpki = 1 ./ (Z .- pol[k])
            hkk = H[k+1, k]
            Q[:, k+1] = (Q[:, k] .* Zpki .- Q[:, 1:k] * H[1:k, k]) ./ hkk
            D[:, k+1] = (D[:, k] .* Zpki .- D[:, 1:k] * H[1:k, k]
                         - Q[:, k] .* Zpki .^ 2) ./ hkk
        end
        R0 = [R0 Q[:, 2:end]]
        R1 = [R1 D[:, 2:end]]
    end
    return R0, R1
end




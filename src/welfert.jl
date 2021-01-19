using Printf

"""
initializes x, d, dp arrays
    input:            n problem size
                      p order of differentiation
                      x[n] collocation points

    output: clog[n]   log of diagonal scaling
            s[n]      sign of diagonal scaling
            d[p*n]    initialized digonal coefficients of Dm

    other : a[n]      alpha(x[n])
            dp[p*n*n] work array

    log of function alpha
"""
function init!(n, p, x, a, clog, s, d, dp)
    
    for i = 1:n
        a[i] = alpha(0, x[i])
    end

    # scaling matrix
    for i = 1:n
        ci = 0.0e0
        si = 1.0e0
        xi = x[i]
        for k = 1:n
            (k == i) && continue
            xik = xi - x[k]
            if xik < 0.0e0
                si = - si
            end
            ci = ci + log(abs(xik))
        end
        clog[i] = ci + a[i]
        s[i] = si
    end

    # initialize d

    for k = 1:n
        xk = x[k]
        for m = 1:p
            i = k + (m - 1) * n
            d[i] = alpha(m, xk)
        end
    end

    return
end



"""
computes diagonal elements of Dm, m = 1[1]p

    input : n      problem size
            p      order of differentiation
            x[n]   collocation points

    output: d[p*n] diagonal elements of Dm
"""
function diag!(n, p, x, d)
    for k = 1:n
        xk = x[k]
        j = k + 1 > n ? 1 : k + 1
        while j != k
            xkj = 1.0e0 / (xk - x[j])
            for m = p:-1:2
                k1 = k + (m - 1) * n
                d[k1] = d[k1] + Float64(m) * xkj * d[k1 - n]
            end
            d[k] = d[k] + xkj
            j = mod(j, n) + 1
        end
    end
    return
end


"""
computes (off-diagonal) elements of Dm, m = 1[1]p

    input : n          problem size
            p          order of differentiation
            x[n]       collocation points
            d[p*n]     diagonal elements of Dm
                       output from subroutine diag

    output: dp[p*n*n]  differentiation matrices Dm

"""
function offdiag!(n, p, x, d, dp)
    for j = 1:n
        j1 = j - 1
        xj = x[j]
        for k = 1:n
            (k == j) && continue 
            i = k + j1 * n
            dp[i] = 1.0e0 / (x[k] - xj)
        end
        dp[j + j1 * n] = d[j]
    end

    (p == 1) && return

    n2 = n * n
    for m = 2:p
        m1 = (m - 2) * n
        m2 = m1 + n
        fm = Float64(m)
        for j = 1:n
            j1 = (m1 + j - 1) * n
            j2 = j1 + n2
            xj = x[j]
            for k = 1:n
                (k == j) && continue 
                k1 = j1 + k
                k2 = j2 + k
                dp[k2] = fm * (d[k + m1] - dp[k1]) / (x[k] - xj)
            end
            dp[j + j2] = d[j + m2]
        end
    end
    return
end


"""
Given n, p, (x_i, i = 1,n),
      and alpha^(j)(x_i), i = 1,n, j = 0,p-1

  computes D1_m and C such that D_m = C D1_m C , m = 1[1]p
  dp[m-1]*n*n+1,...,m*n*n) contains D1_m
  Cii = s[i]*exp(clog[i]), i = 1[1]n
"""
function  make_scaled_dp()
    #integer nmax, pmax
    nmax = 64
    pmax = 8
    
    x    = zeros(nmax)
    d    = zeros(nmax * pmax)
    dp   = zeros(pmax * nmax * nmax)
    a    = zeros(nmax)
    clog = zeros(nmax)
    s    = zeros(nmax)

    p = 2
    n = 8
    
    # collocation points
    for i ∈ 1:n
        x[i] = 10.0e0 * (i-1) / (n-1)
    end

    init!(n, p, x, a, clog, s, d, dp)
    diag!(n, p, x, d)
    offdiag!(n, p, x, d, dp)

    istart = (p - 1) * n * n + 1
    iend = p * n * n

    # to match printing from the Welfert paper
    for i ∈ istart:8:iend
        @printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", dp[i:i+7]...)
    end

    @printf("\n%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", s[1:8]...)
    @printf("\n%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", clog[1:8]...)

end    


"""
    if m = 0: log of alpha(x)
    if m > 0: alpha^m(x)/alpha(x)

    Here alpha(x) = 1 / (1 + x) ** n
"""
function alpha(m, x)
    n = 8
    if m == 0
        alpha = - float(n) * log(1.0e0 + x)
    else
        ai = 1.0e0
        for i = n:n + m - 1
            ai = - ai * float(i) / (1.0e0 + x)
        end
        alpha = ai
    end
    return alpha
end

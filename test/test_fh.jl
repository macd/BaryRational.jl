
runge(x) = 1 / (1 + 25 * x^2)
    
function test_fh_runge(tol=1e-8)
    xx = [-1.0:0.05:1.0;]
    yy = runge.(xx)
    fh = FHInterp(xx, yy, order=3)
    xxx = [-1.0:0.01:1.0;]
    vl = fh.(xxx)
    rf = runge.(xxx)
    return norm(vl - rf) < 5e-5
end

# Ok, so this is kinda brittle..
function test_abs_x(order=3)
    pts = [10, 20, 40, 80, 160, 320, 640]
    tol = [0.5, 0.2, 0.09, 0.05, 0.03, 0.01, 0.005]
    pass = true
    for (p, t) in zip(pts, tol)
        xx = collect(range(-5.0, 5.0, length=2p-1))
        xi = xx[1:2:end]
        xt = xx[2:2:end]
        yy = abs.(xi)
        fh = FHInterp(xi, yy, order=order, grid=true)
        err = maximum(abs.(fh.(xt) .- abs.(xt)))
        pass *= err < t
    end
    return pass
end

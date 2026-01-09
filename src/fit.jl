function doublet_gauss((x, y), (bg, A1, A2, x1, y1, x2, y2, dx, dy))
    bg +
    A1 * exp(-((x-x1)/dx)^2/2 - ((y-y1)/dy)^2/2) +
    A2 * exp(-((x-x2)/dx)^2/2 - ((y-y2)/dy)^2/2)
end

function fit_peak(
    img::AbstractMatrix,
    (x1, y1),
    (x2, y2);
    padding = 10,
    sigma = 1,
    reltol = 1e-15,
)
    xc, yc = round.(Int, [x1 + x2, y1 + y2] / 2)
    box = CartesianIndices((xc-padding:xc+padding, yc-padding:yc+padding))
    boxed = Float64.(img[box])
    bg = median(img)
    Am = maximum(boxed)
    x0, y0 = first(box).I .- 1
    u0 = Float64[bg, Am, Am/2, x1 - x0, y1 - y0, x2 - x0, y2 - y0, sigma, sigma]
    # Am, xym = findmax(boxed)
    # As, xys = if x1 < x2
    #     findmax(boxed[round(Int, xym[1]+3sigma):end, :])
    # else 
    #     findmax(boxed[begin:round(Int, xym[1]-3sigma), :])
    # end
    # u0 = Float64[bg, Am, As, xym.I..., xys.I..., sigma, sigma]
    resid = (u, p) -> sum(CartesianIndices(p)) do i
        abs2(p[i] - doublet_gauss(i.I, u))
    end
    optf = OptimizationFunction(resid, SecondOrder(AutoForwardDiff(), AutoForwardDiff()))
    prob = OptimizationProblem(optf, u0, boxed)
    sol = solve(prob, Newton(alphaguess = 0.001); reltol)
    bg, A1, A2, x1, y1, x2, y2, _, _ = sol.u
    (A1 = A1, x1 = x0 + x1, y1 =y0 + y1, A2 = A2, x2 = x0 + x2, y2 = y0 + y2)
end

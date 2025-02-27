function profile(model::Model, peak::Peak; mul = 1)
    hkl = reflex(model, peak.coord)
    grid = Tuple((range(first(ax), last(ax), length(ax) * mul) for ax in peak.indices))
    vals = [idx -> intensity(model, hkl, idx.I) for idx in product(grid...)]
    linear_interpolation(grid, vals, extrapolation_bc = 0.0)
end

struct FitParams{N}
    shift::Vec{N}
    scale::Number
    sigma::Number
    baseline::Number
end

function profile_fit(coords::AbstractArray, profile, p::FitParams)
    shifted = map(coord -> profile((coord .- p.shift)...), coords)
    p.scale * imfilter(shifted, Kernel.gaussian(p.sigma)) .+ p.baseline
end

function refine_profile(peak::Peak{N}, profile; maxscale = 4, maxsigma = 16) where {N}
    coords = Tuple.(collect(CartesianIndices(peak.image)))
    profile_int = sum(coord -> profile(coord...), coords)
    scale0 = peak.int / profile_int
    baseline0 = median(peak.image)

    lsq = OptimizationFunction(
        (u, p) -> begin
            shift..., scale, sigma, baseline = u
            params = FitParams(shift, scale, sigma, baseline)
            vals = profile_fit(coords, profile, params)
            sum(abs2, image .- vals)
        end,
        AutoFiniteDiff(),
    )
    solve(
        OptimizationProblem(
            lsq,
            [zeros(N)..., scale0, 1.0, baseline0],
            lb = [-ones(N)..., scale0 / maxscale, 0.0, minimum(peak.image)],
            ub = [onces(N)..., scale0 * maxscale, maxsigma, mean(peak.image)],
        ),
        Optimization.LBFGS(),
    ).u
end

function refined_position(peak::Peak{N}, params::FitParams) where {N}
    Vec{N}(peak.coord .+ params.shift)
end

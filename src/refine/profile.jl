function profile_correction(coords::CartesianIndices, profile, p)
    shift..., scale, sigma, baseline = p
    shifted = map(coord -> profile((coord .- shift)...), coords)
    exp(scale) * imfilter(shifted, Kernel.gaussian(sigma)) .+ baseline
end

function refine_peak(
    experiment::Experiment,
    peak::Peak{N};
    mul = 2,
    maxshift = 3,
    maxscale = 10,
    maxsigma = 16,
) where {N}
    model = model(experiment, peak.frame.setting)
    profile = profile(model, peak.region.indices, peak.coord; mul = mul)
    func = p -> profile_correction(peak.region, profile, p)
    image = view(peak.frame.image, peak.region)

    p0 = [zeros(N)..., 1.0, 1.0, 0.0]
    profile_intensity = sum(func, p0)
    scale0 = peak.intensity / profile_intensity
    p0 = [zeros(N)..., scale0, 1.0, peak.baseline]
    pmin = [fill(maxshift, N)..., scale0 * maxscale, 0.0, 0.0]
    pmax = [fill(-maxshift, N)..., scale0 / maxscale, maxsigma, 3 * peak.baseline]

    lsq = OptimizationFunction((p, _) -> sum(abs2, image .- func(p)), AutoFiniteDiff())
    fit = solve(OptimizationProblem(lsq, p0, lb = pmin, ub = pmax), Optimization.LBFGS()).u

    shift..., scale, sigma, baseline = fit
    Peak(
        sum(func, [shift..., scale, sigma]),
        baseline,
        peak.coord .+ shift,
        peak.region,
        peak.frame,
    )
end

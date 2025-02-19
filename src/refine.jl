function refine_profile(
    image::AbstractArray,
    coords::AbstractArray,
    profile::Profile{N},
) where {N}
    image_stat = sum((coord, int) -> [1, coord...] * int, zip(coords, image))
    scale = image_stat[1] / intensity(profile)
    shift = image_stat[2:end] / image_stat[1] .- coord_mean(profile)
    trans = Mat{N,N}(I)
    lsq = OptimizationFunction(
        (params, _) -> mean(abs2, image .- AffinedProfile(profile, params...).(coords)),
    )
    fit_params = solve(OptimizationProblem(lsq, [scale, shift..., trans...])).u
    AffinedProfile(profile, fit_params...)
end

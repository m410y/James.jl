function profile_correction(
    scale::Number,
    shift_x::Number,
    shift_y::Number,
    sigma_x::Number,
    sigma_y::Number,
    baseline::Number,
    coords::AbstractArray,
    profile,
)
    shifted = map(coord -> profile(coord[1] - shift_x, coord[2] - shift_y), coords)
    scale * imfilter(shifted, Kernel.gaussian((sigma_x, sigma_y))) .+ baseline
end

function refine_profile(image::AbstractArray, coords::AbstractArray, profile)
    lsq = OptimizationFunction(
        (u, p) -> begin
            vals = profile_correction(u..., coords, profile)
            sum(abs2, image .- vals)
        end,
        AutoFiniteDiff(),
    )
    solve(
        OptimizationProblem(
            lsq,
            [0.5, 0.0, -13.0, 1.0, 1.0, 80.0],
            lb = [0.1, -10.0, -20.0, 0.0, 0.0, 0.0],
            ub = [10.0, 10.0, 0.0, 3.0, 3.0, 200.0],
        ),
        Optimization.LBFGS(),
    ).u
end

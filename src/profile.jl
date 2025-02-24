function profile(model::Model, hkl::AbstractVector, grid; baseline = 0)
    vals = [intensity(model, hkl, coord) for coord in product(grid...)]
    itp = interpolate(grid, vals, Gridded(Linear()))
    extrapolate(itp, baseline)
end

struct SimpleConvParams{N}
    scale::Number
    shift::Vec{N}
    sigma::Vec{N}
    baseline::Number
end

function profile_convolution(
    params::SimpleConvParams,
    coords::AbstractArray,
    profile,
)
    shifted = map(coord -> profile(coord[1] - shift_x, coord[2] - shift_y), coords)
    scale * imfilter(shifted, Kernel.gaussian((sigma_x, sigma_y))) .+ baseline
end

function image_stat(image::AbstractArray, coords::AbstractArray)
    stat = sum(p -> p[1] * [1, p[2]...] * [1, p[2]...]', zip(image, coords))
    stat0 = stat[begin]
    stat1 = stat[begin+1:end, begin]
    stat2 = stat[begin+1:end, begin+1:end]
    stat0, stat1 / stat0, (stat2 * stat0 - stat1 * stat1') / stat0^2
end

function refine_profile(image::AbstractArray, coords::AbstractArray, profile)
    imstat = image_stat(image, coords)

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



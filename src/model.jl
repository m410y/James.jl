abstract type Model end

struct StillModel{K<:Spectrum,D<:Detector,S<:Sample} <: Model
    spectrum::K
    detector::D
    sample::S
end

struct ScanModel{K<:Spectrum,D<:Detector,S<:Sample} <: Model
    spectrum::K
    detector::D
    sample::S
    axis::Axis
    increment::Number
end

StillModel(model::ScanModel, angle::Number) =
    StillModel(model.spectrum, model.detector, model.axis(angle)(model.sample))

function reflection_angle(model::ScanModel, hkl::AbstractVector)
    k0 = mean(model.spectrum)
    s0 = model.sample.UB * hkl
    angles = reflection_angles(model.axis, s0, k0)
    nearest = findmin(abs, rem2pi.(NoUnits.(angles .- model.increment / 2), RoundNearest))
    angles[last(nearest)]
end

function intensity(model::StillModel, hkl::AbstractVector, coord::AbstractVector)
    s = model.sample.UB * hkl
    r = model.detector(coord) - model.sample.p
    n = normalize(r)
    k = n * dot(s, s) / 2dot(n, s) - s
    model.spectrum(k)
end

function intensity(model::ScanModel, hkl::AbstractVector, coord::AbstractVector)
    @show angle = reflection_angle(model, hkl)
    @show sample0 = model.axis(angle)(model.sample)
    @show s = sample0.UB * hkl
    @show r = model.detector(coord) - sample0.p
    @show n = normalize(r)
    @show k = n * dot(s, s) / 2dot(n, s) - s
    @show ds = cross(model.axis.v, s)
    @show dk = n * dot(s, s) * dot(n, ds) / 2dot(n, s)^2 - ds

    @show A = inv(cov(model.spectrum))
    @show c1 = dot(dk, A, dk)
    @show c2 = dot(dk, A, k - mean(model.spectrum))
    @show C = Symmetric(Mat2(c1, c2, c2, 1))
    @show eig = eigen(C)
    if det(eig) > 0
        return 0
    end
    N = eig.vectors * [0 1; 1 0] * sqrt(abs.(Diagonal(eig.values))) * [1 1; 1 -1]
    lb, ub = extrema(N[1, :] ./ N[2, :])
    if lb > model.increment || ub < 0
        return 0
    end

    return 1
    problem = IntegralProblem((u, _) -> model.spectrum(k + dk * u), (lb, ub))
    solve(problem, HCubatureJL()).u
end

function reflex(model::StillModel, coord::AbstractVector)
    r = model.detector(coord...) - model.sample.p
    n = normalize(r)
    k0 = mean(model.spectrum)
    s = norm(k0) * n - k0
    model.sample.UB \ s
end

reflex(model::ScanModel, coord::AbstractVector) =
    reflex(StillModel(model, model.increment / 2), coord)

function coord(model::StillModel, hkl::AbstractVector)
    k0 = mean(model.spectrum)
    s = model.sample.UB * hkl
    kd = k0 + s * (1 / 2 - dot(k0, s) / dot(s, s))
    intersect_coord(model.detector, model.sample.p, kd)
end

function coord(model::ScanModel, hkl::AbstractVector)
    k0 = mean(model.spectrum)
    s0 = model.sample.UB * hkl
    angle = reflection_angle(model, hkl)
    kd = k0 + model.axis(angle)(s0)
    intersect_coord(model.detector, model.sample.p, kd)
end

function profile(model::Model, indices::NTuple{N,AbstractRange}, coord; mul = 1) where {N}
    hkl = reflex(model, coord)
    grid = Tuple((range(first(ax), last(ax), length(ax) * mul) for ax in indices))
    vals = [intensity(model, hkl, idx) for idx in product(grid...)]
    linear_interpolation(grid, vals, extrapolation_bc = 0.0)
end

function Base.show(io::IO, ::MIME"text/plain", model::StillModel)
    show(io, "text/plain", model.spectrum)
    println(io)
    show(io, "text/plain", model.detector)
    println(io)
    show(io, "text/plain", model.sample)
end

function Base.show(io::IO, ::MIME"text/plain", model::ScanModel)
    show(io, "text/plain", model.spectrum)
    println(io)
    show(io, "text/plain", model.detector)
    println(io)
    show(io, "text/plain", model.sample)
    println(io)
    show(io, "text/plain", model.axis)
    println(io)
    print(io, "  increment: $(model.increment)")
end

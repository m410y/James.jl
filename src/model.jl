abstract type Model end

struct StillModel{K<:Spectrum,D<:Detector,S<:Sample} <: Model
    spectrum::K
    detector::D
    sample::S
end

struct ScanModel{T,K<:Spectrum,D<:Detector,S<:Sample} <: Model
    spectrum::K
    detector::D
    sample::S
    axis::Axis{T}
    increment::T
end

StillModel(model::ScanModel, angle) =
    StillModel(model.spectrum, model.detector, model.axis(angle)(model.sample))

function intensity(model::StillModel, hkl, coord)
    s = model.sample.UB * Vec(hkl...)
    r = model.detector(coord...) - model.sample.p
    n = normalize(r)
    k = n * dot(s, s) / 2dot(n, s) - s
    model.spectrum(k)
end

function intensity(model::ScanModel, hkl, coord)
    problem = IntegralProblem(
        (angle, _) -> intensity(StillModel(model, angle), hkl, coord),
        (0, model.increment),
    )
    solve(problem, HCubatureJL()).u
end

function reflex(model::StillModel, coord)
    r = model.detector(coord...) - model.sample.p
    n = normalize(r)
    k0 = mean(model.spectrum)
    s = norm(k0) * n - k0
    inv(Matrix(model.sample.UB)) * s
end

reflex(model::ScanModel, coord) = reflex(StillModel(model, model.increment / 2), coord)

function coord(model::StillModel, hkl)
    k0 = mean(model.spectrum)
    s = model.sample.UB * Vec(hkl...)
    kd = k0 + s * (1 / 2 - dot(k0, s) / dot(s, s))
    intersect_coord(model.detector, model.sample.p, kd)
end

function coord(model::ScanModel, hkl)
    k0 = mean(model.spectrum)
    s0 = model.sample.UB * Vec(hkl...)
    angles = reflection_angles(model.axis, s0, k0)
    nearest = findmin(abs, rem2pi.(NoUnits.(angles .- model.increment / 2), RoundNearest))
    angle = angles[last(nearest)]
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

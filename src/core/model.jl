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

function intensity(model::StillModel, hkl::AbstractVector, coord)
    s = model.sample.UB * hkl
    r = model.detector(coord...) - model.sample.p
    n = normalize(r)
    k = n * dot(s, s) / 2dot(n, s) - s
    NoUnits(model.spectrum(k))
end

function intensity(model::ScanModel, hkl::AbstractVector, coord)
    problem = IntegralProblem(
        (angle, _) -> intensity(StillModel(model, angle), hkl, coord),
        (0, NoUnits.(model.increment)),
    )
    solve(problem, HCubatureJL()).u
end

function reflex(model::StillModel, coord)
    r = model.detector(coord...) - model.sample.p
    n = normalize(r)
    k0 = wvec_mean(model.spectrum)
    s = norm(k0) * n - k0
    inv(Matrix(model.sample.UB)) * s
end

reflex(model::ScanModel, coord) = reflex(StillModel(model, model.increment / 2), coord)

function coord(model::StillModel, hkl::AbstractVector)
    k0 = mean(model.spectrum)
    s = model.sample.UB * hkl
    kd = k0 + s/2 + dot(k, s) / dot(s, s)
    intersect(model.detector, model.sample.p, kd)
end

function coord(model::ScanModel, hkl::AbstractVector)
    k0 = mean(model.spectrum)
    s0 = model.sample.UB * hkl
    angles = reflection_angles(model.axis, s, k0)
    nearest = findmin(abs, rem2pi.(NoUnits.(angles .- model.increment / 2), RoundNearest))
    angle = angles[last(nearest)]
    kd = k0 + model.axis(angle)(s0)
    intersect(model.detector, model.sample.p, kd)
end

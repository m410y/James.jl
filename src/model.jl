abstract type Model{N} end

struct StillModel{N} <: Model{N}
    spec::Spectrum
    detect::Detector{N}
    cryst::SingleCrystal
end

struct ScanModel{N} <: Model{N}
    spec::Spectrum
    detect::Detector{N}
    cryst::SingleCrystal
    axis::Axis
    angles::NTuple{2,Number}
end

StillModel(model::ScanModel, angle::Number) =
    StillModel(model.spec, model.detect, model.axis(angle)(model.cryst))

function intensity(model::StillModel, hkl::AbstractVector, coord)
    s = model.cryst.UB * hkl
    r = model.detect(coord...) - model.cryst.pos
    n = normalize(r)
    k = n * dot(s, s) / 2dot(n, s) - s
    NoUnits(model.spec(k))
end

function intensity(model::ScanModel, hkl::AbstractVector, coord)
    problem = IntegralProblem(
        (angle, _) -> intensity(StillModel(model, angle), hkl, coord),
        NoUnits.(model.angles)
    )
    solve(problem, HCubatureJL()).u
end

function profile(model::StillModel{N}, hkl::AbstractVector, grid; baseline = 0) where {N}
    vals = map(coord -> intensity(model, hkl, coord), product(grid...))
    itp = interpolate(grid, vals, Gridded(Linear()))
    GriddedProfile{N}(extrapolate(itp, baseline))
end

function profile(model::ScanModel{N}, hkl::AbstractVector, grid; baseline = 0) where {N}
    vals = [intensity(model, hkl, coord) for coord in product(grid...)]
    itp = interpolate(grid, vals, Gridded(Linear()))
    GriddedProfile{N}(extrapolate(itp, baseline))
end

function reflex(model::StillModel, coord)
    r = model.detect(coord...) - model.cryst.pos
    n = normalize(r)
    k0 = wvec_mean(model.spec)
    s = norm(k0) * n - k0
    Vec3(inv(Matrix(model.cryst.UB)) * s)
end

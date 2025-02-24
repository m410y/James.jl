abstract type Model{N} end

struct MovingModel{K,M,N} <: Model{N}
    spec::Spectrum
    detect::Detector{N}
    gonio_d::Goniometer{M}
    cryst::SingleCrystal
    gonio_c::Goniometer{K}
end

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

StillModel(model::MovingModel, angles_d, angles_c) = 
    StillModel(
        model.spec,
        model.gonio_d(angles_d...)(model.detect),
        model.gonio_c(angles_c...)(model.cryst),
    )

ScanModel(model::MovingModel, angles_d, angles_c, (n_axis, inc)) = 
    ScanModel(
        model.spec,
        model.gonio_d(angles_d...)(model.detect),
        model.gonio_c(angles_c...)(model.cryst),
        model.gonio_c.axes[n_axis],
        (0.0, inc)
    )

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
        NoUnits.(model.angles),
    )
    solve(problem, HCubatureJL()).u
end

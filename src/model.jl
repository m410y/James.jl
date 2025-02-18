function intensity_still_precise(
    spec::Spectrum,
    detector::Detector2D,
    cryst::SingleCrystal,
    hkl::AbstractVector,
    x::Number,
    y::Number,
)
    s = reflex(cryst, hkl)
    n = Vec3(normalize(detector(x, y) - cryst.pos))
    k = n * dot(s, s) / 2dot(n, s) - s
    spec(k)
end

function profile_still_approx(
    spec::Spectrum,
    detector::Detector2D,
    cryst::SingleCrystal,
    hkl::AbstractVector,
    x::Number,
    y::Number,
)
    s = reflex(cryst, hkl)
    r = detector(x, y) - cryst.pos
    n, d = Vec3(normalize(r)), norm(r)
    k = n * dot(s, s) / 2dot(n, s) - s
    nx = (detector.ex - n * dot(n, detector.ex)) / d
    ny = (detector.ey - n * dot(n, detector.ey)) / d
    kx = dot(s, s) / 2dot(n, s) * (nx - n * dot(s, nx) / dot(n, s))
    ky = dot(s, s) / 2dot(n, s) * (ny - n * dot(s, ny) / dot(n, s))
    SpectrumSlice(spec, Vec3(k - kx * x - ky * y), Vec3(kx), Vec3(ky))
end

function profile_scan_approx(
    spec::Spectrum,
    detector::Detector2D,
    (cryst_start, cryst_stop)::Pair{SingleCrystal,SingleCrystal},
    hkl::AbstractVector,
    x::Number,
    y::Number,
)
    s_start = reflex(cryst_start, hkl)
    s_stop = reflex(cryst_stop, hkl)
    s = (s_start + s_stop) / 2
    ds = (s_stop - s_start) / 2
    r = detector(x, y) - (cryst_start.pos + cryst_stop.pos) / 2
    n, d = normalize(r), norm(r)
    k = n * dot(s, s) / 2dot(n, s) - s
    nx = (detector.ex - n * dot(n, detector.ex)) / d
    ny = (detector.ey - n * dot(n, detector.ey)) / d
    kx = dot(s, s) / 2dot(n, s) * (nx - n * dot(s, nx) / dot(n, s))
    ky = dot(s, s) / 2dot(n, s) * (ny - n * dot(s, ny) / dot(n, s))
    kv = n * (dot(s, ds) / dot(n, s) - dot(s, s) * dot(n, ds) / 2dot(n, s)^2) - ds
    SpectrumProjection(spec, Vec3(k - kx * x - ky * y), Vec3(kx), Vec3(ky), Vec3(kv))
end

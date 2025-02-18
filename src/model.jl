function intensity_still(
    spec::Spectrum,
    detect::Detector,
    cryst::SingleCrystal,
    hkl::AbstractVector,
    coord,
)
    s = reflex(cryst, hkl)
    r = Vec3(detect(coord...) - cryst.pos)
    n, d = normalize(r), norm(d)
    k = n * dot(s, s) / 2dot(n, s) - s
    spec(k) / d^2
end

function intensity_scan(
    spec::Spectrum,
    detect::Detector,
    cryst::SingleCrystal,
    hkl::AbstractVector,
    coord,
    axis::Axis,
    angles::Pair{Number,Number},
)
    solve(
        IntegralProblem(
            (u, p) -> intensity_still(spec, detect, axis(u)(cryst), hkl, coord),
            angles,
        ),
        HCubatureJL()
    ).u
end

function intensity_still(
    spec::Spectrum,
    detect::Detector,
    cryst::SingleCrystal,
    hkl::AbstractVector,
    coord,
)
    s = reflex(cryst, hkl)
    r = Vec3(detect(coord...) - cryst.pos)
    n, d = normalize(r), norm(r)
    k = n * dot(s, s) / 2dot(n, s) - s
    spec(k)
end

function intensity_scan(
    spec::Spectrum,
    detect::Detector,
    cryst::SingleCrystal,
    hkl::AbstractVector,
    coord,
    axis::Axis,
    angles::NTuple{2,Number},
)
    solve(
        IntegralProblem(
            (u, p) -> intensity_still(spec, detect, axis(u)(cryst), hkl, coord),
            NoUnits.(angles),
        ),
        HCubatureJL()
    ).u
end

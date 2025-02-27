abstract type Spectrum end

struct DoubletKαX <: Spectrum
    I::Number
    ratio::Number
    E1::Energy
    E2::Energy
    ΔE::Energy
    Δk::Energy
end

mean(spec::DoubletKαX) =
    Vec3([1, 0, 0]) * (spec.ratio * spec.E1 + spec.E2) / (spec.ratio + 1)

function (spec::DoubletKαX)(k)
    divergency = exp(-(k[2]^2 + k[3]^2) / spec.Δk^2 / 2) / (2π * spec.Δk^2)
    peak1 = lorentz(k[1], spec.E1, spec.ΔE)
    peak2 = lorentz(k[1], spec.E2, spec.ΔE)
    peak_sum = (spec.ratio * peak1 + peak2) / (spec.ratio + 1)
    spec.I * divergency * peak_sum
end

struct GriddedSpectrum <: Spectrum
    k0::Vec3
    M::Mat3
    itp::Interpolations.Extrapolation
end

mean(spec::GriddedSpectrum) = spec.k0

function (spec::GriddedSpectrum)(k)
    itp_coords = M * (k - k0)
    itp(itp_coords...)
end

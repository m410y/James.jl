abstract type Spectrum end

function intensity(spec::Spectrum)
    solve(
        IntegralProblem((k, p) -> spec(k), [-Inf, -Inf, -Inf], [Inf, Inf, Inf]),
        HCubatureJL(),
    ).u
end

function wvec_mean(spec::Spectrum)
    int =
        solve(
            IntegralProblem(
                (k, p) -> [k; 1] * spec(k),
                [-Inf, -Inf, -Inf],
                [Inf, Inf, Inf],
            ),
            HCubatureJL(),
        ).u
    int[1:3] / int[4]
end

function wvec_disp(spec::Spectrum)
    int =
        solve(
            IntegralProblem(
                (k, p) -> [k; 1] * [k; 1]' * spec(k),
                [-Inf, -Inf, -Inf],
                [Inf, Inf, Inf],
            ),
            HCubatureJL(),
        ).u
    (int[1:3, 1:3] * int[4, 4] - int[1:3, 4] * int[4, 1:3]') / int[4, 4]^2
end

struct DoubletKαX <: Spectrum
    I::Number
    ratio::Number
    E1::Number
    E2::Number
    ΔE::Number
    Δk::Number
end

function (spec::DoubletKαX)(k::AbstractVector)
    @assert length(k) == 3
    divergency = exp(-(k[2]^2 + k[3]^2) / spec.Δk^2 / 2) / (2π * spec.Δk^2)
    peak1 = lorentz(k[1], spec.E1, spec.ΔE)
    peak2 = lorentz(k[1], spec.E2, spec.ΔE)
    peak_sum = (spec.ratio * peak1 + peak2) / (spec.ratio + 1)
    spec.I * divergency * peak_sum
end

intensity(spec::DoubletKαX) = spec.I
wvec_mean(spec::DoubletKαX) =
    Vec3([1, 0, 0]) * (spec.ratio * spec.E1 + spec.E2) / (spec.ratio + 1)
wvec_disp(spec::DoubletKαX) = Diagonal([Inf * spec.ΔE^2, spec.Δk^2, spec.Δk^2])

struct GriddedSpectrum
    ext::Interpolations.Extrapolation
end

(grid::GriddedSpectrum)(k::AbstractVector) = ext(k)

abstract type Spectrum{T} end

struct DoubletKαX{T} <: Spectrum{T}
    I::T
    ratio::T
    E1::T
    E2::T
    ΔE::T
    Δk::T
end

mean(spec::DoubletKαX) =
    Vec3([1, 0, 0]) * (spec.ratio * spec.E1 + spec.E2) / (spec.ratio + 1)

function _gauss(x::Number, x0::Number, σ::Number)
    exp(-((x - x0) / σ)^2 / 2) / sqrt(2π) / σ
end

function _lorentz(x::Number, x0::Number, γ::Number)
    1 / (1 + ((x - x0) / γ)^2) / π / γ
end

function (spec::DoubletKαX)(k)
    divergency = exp(-(k[2]^2 + k[3]^2) / spec.Δk^2 / 2) / (2π * spec.Δk^2)
    peak1 = _lorentz(k[1], spec.E1, spec.ΔE)
    peak2 = _lorentz(k[1], spec.E2, spec.ΔE)
    peak_sum = (spec.ratio * peak1 + peak2) / (spec.ratio + 1)
    spec.I * divergency * peak_sum
end

function Base.show(io::IO, ::MIME"text/plain", spec::DoubletKαX)
    print(io, "Simple Kα spectrum along X axis:\n")
    print(io, "  intensity: $(spec.I)\n")
    print(io, "  intensity ratio: $(spec.ratio)\n")
    print(io, "  Kα₁ energy: $(spec.E1 * ReciprocalUnit)\n")
    print(io, "  Kα₂ energy: $(spec.E2 * ReciprocalUni)\n")
    print(io, "  X axis width: $(spec.ΔE * ReciprocalUni)\n")
    print(io, "  YZ plane width: $(spec.Δk * ReciprocalUni)")
end


struct GriddedSpectrum{T} <: Spectrum{T}
    k0::Vec3{T}
    M::Mat3{T}
    itp::Interpolations.Extrapolation{T,3}
end

mean(spec::GriddedSpectrum) = spec.k0

function (spec::GriddedSpectrum)(k)
    itp_coords = M * (Vec(k...) - k0)
    itp(itp_coords...)
end

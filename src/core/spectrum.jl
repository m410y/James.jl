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

function (spec::DoubletKαX)(k::AbstractVector)
    divergency = exp(-(k[2]^2 + k[3]^2) / spec.Δk^2 / 2) / (2π * spec.Δk^2)
    peak1 = _lorentz(k[1], spec.E1, spec.ΔE)
    peak2 = _lorentz(k[1], spec.E2, spec.ΔE)
    peak_sum = (spec.ratio * peak1 + peak2) / (spec.ratio + 1)
    spec.I * divergency * peak_sum
end
(spec::DoubletKαX)(k::Vararg{Number,3}) = spec(Vec3(k))

function Base.show(io::IO, ::MIME"text/plain", spec::DoubletKαX; unit = u"eV")
    convert(E) = unit != ReciprocalUnit ? ustrip(uconvert(unit, E * ReciprocalUnit, Spectral())) : E
    E0 = norm(mean(spec))
    println(io, summary(spec), ":")
    println(io, "  intensity: $(spec.I)")
    println(io, "  intensity ratio: $(spec.ratio)")
    println(io, "  Kα₁ energy: ", @sprintf("%8.3f", convert(spec.E1)), " $unit")
    println(io, "  Kα₂ energy: ", @sprintf("%8.3f", convert(spec.E2)), " $unit")
    println(io, "  X width   : ", @sprintf("%8.3f", abs(convert(E0 + spec.ΔE) - convert(E0))), " $unit")
    println(io, "  YZ width  : ", @sprintf("%8.3f", abs(convert(E0 + spec.Δk) - convert(E0))), " $unit")
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

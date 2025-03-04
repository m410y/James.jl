abstract type Spectrum{T} end

struct DistSpectrum{T,D<:Distribution} <: Spectrum{T}
    intensity::T
    mean::Vec3{T}
    cov::Symmetric{T,Mat3{T}}
    dist::D
    function DistSpectrum(
        intensity::Number,
        mean::AbstractVector,
        cov::AbstractMatrix,
        dist::D,
    ) where {D<:Distribution}
        new{Float64,D}(intensity, Vec3d(mean), Symmetric(Mat3d(cov)), dist)
    end
end

function DistSpectrum(
    intensity::Number,
    dist::Distribution,
    domain::NTuple{2,AbstractVector};
    alg = HCubatureJL(),
)
    ifunc(u, _) = [u..., 1] * [u..., 1]' * intensity * pdf(dist, u)
    solution = solve(IntegralProblem(ifunc, domain), alg)
    if solution.retcode != ReturnCode.Success
        error("Cant integrate function")
    end
    intensity = solution.u[4, 4]
    mean = solution.u[1:3, 4] / intensity
    cov = solution.u[1:3, 1:3] / intensity - mean * mean'
    DistSpectrum(intensity, mean, cov, fun)
end

(spec::DistSpectrum)(k::AbstractVector) = spec.intensity * pdf(spec.dist, k)
intensity(spec::DistSpectrum) = spec.intensity
Statistics.mean(spec::DistSpectrum) = spec.mean
Statistics.cov(spec::DistSpectrum) = spec.cov

struct SpectrumSum{T,N} <: Spectrum{T}
    intensity::T
    mean::Vec3{T}
    cov::Symmetric{T,Mat3{T}}
    specs::NTuple{N,Spectrum}
    function SpectrumSum(
        intensity::Number,
        mean::AbstractVector,
        cov::AbstractMatrix,
        dist::NTuple{N,Spectrum},
    ) where {N}
        new{Float64,N}(intensity, Vec3d(mean), Symmetric(Mat3d(cov)), dist)
    end
end

function SpectrumSum(specs::Vararg{Spectrum,N}) where {N}
    Q = sum(s -> let I = intensity(s), E = mean(s), C = cov(s)
        [C+E*E' E; E' 1] * I
    end, specs)
    sum_intensity = Q[4, 4]
    sum_mean = Q[1:3, 4] / sum_intensity
    sum_cov = Q[1:3, 1:3] / sum_intensity - sum_mean * sum_mean'
    SpectrumSum(sum_intensity, sum_mean, sum_cov, specs)
end

(spec::SpectrumSum)(k::AbstractVector) = sum(s -> s(k), spec.specs)
intensity(spec::SpectrumSum) = spec.intensity
Statistics.mean(spec::SpectrumSum) = spec.mean
Statistics.cov(spec::SpectrumSum) = spec.cov

function Base.show(io::IO, ::MIME"text/plain", spec::Spectrum)
    kunit = u"eV"
    println(summary(spec), ":")
    println("  intensity: ", round(intensity(spec), sigdigits = 5))
    println(
        "  mean [$kunit]: ",
        @sprintf("%8.2f, %8.2f, %8.2f", reciprocal_convert.(kunit, mean(spec))...)
    )
    println("  covariation [$kunit]:")
    for row in eachrow(reciprocal_convert.(kunit, cov(spec)))
        println(io, @sprintf("    %8.2f %8.2f %8.2f", row...))
    end
end

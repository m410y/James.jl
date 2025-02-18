abstract type Profile2D end

function intensity(prof::Profile2D)
    solve(
        IntegralProblem(((x, y), p) -> prof(x, y), [-Inf, -Inf], [Inf, Inf]),
        HCubatureJL(),
    ).u
end

function coord_mean(prof::Profile2D)
    int =
        solve(
            IntegralProblem(
                ((x, y), p) -> [x, y, 1] * prof(x, y),
                [-Inf, -Inf],
                [Inf, Inf],
            ),
            HCubatureJL(),
        ).u
    int[1:2] / int[3]
end

function coord_disp(prof::Profile2D)
    int =
        solve(
            IntegralProblem(
                ((x, y), p) -> [x, y, 1] * [x, y, 1]' * prof(x, y),
                [-Inf, -Inf],
                [Inf, Inf],
            ),
            HCubatureJL(),
        ).u
    (int[1:2, 1:2] * int[3, 3] - int[1:2, 3] * int[3, 1:2]') / int[3, 3]^2
end

struct Gaussian <: Profile2D end
(::Gaussian)(x::Number, y::Number) = exp(-(x^2 + y^2) / 2) / 2π
intensity(::Gaussian) = 1
coord_mean(::Gaussian) = [0, 0]
coord_disp(::Gaussian) = Diagonal([1, 1])

struct SpectrumSlice{S<:Spectrum} <: Profile2D
    spec::S
    k0::Vec3
    kx::Vec3
    ky::Vec3
end

(slice::SpectrumSlice)(x::Number, y::Number) =
    slice.spec(slice.k0 + x * slice.kx + y * slice.ky)

struct SpectrumProjection{S<:Spectrum} <: Profile2D
    spec::S
    k0::Vec3
    kx::Vec3
    ky::Vec3
    kv::Vec3
end

function (proj::SpectrumProjection)(x::Number, y::Number)
    solve(
        IntegralProblem(
            (t, p) -> 1e15*ustrip(proj.spec(proj.k0 + x * proj.kx + y * proj.ky + t * proj.kv)),
            (-1, 1),
        ),
        HCubatureJL(), reltol=1e-10
    ).u
end

struct GriddedProfile <: Profile2D
    ext::Interpolations.AbstractExtrapolation
end

(grid::GriddedProfile)(x::Number, y::Number) = grid.ext(x, y)

struct AffinedProfile{P<:Profile2D} <: Profile2D
    profile::P
    A::Number
    M::Mat2
    V::Vec2
end

function (aff::AffinedProfile)(x::Number, y::Number)
    new_x, new_y = aff.M * [x, y] + aff.V
    aff.A * aff.profile(new_x, new_y)
end
intensity(aff::AffinedProfile) = aff.A / det(aff.M) * intensity(aff.profile)
coord_mean(aff::AffinedProfile) = inv(aff.M) * (coord_mean(aff.profile) - aff.V)
coord_disp(aff::AffinedProfile) = inv(aff.M) * coord_disp(aff.profile) * inv(aff.M)'

struct ProfileSum{N} <: Profile2D
    profs::NTuple{N,Profile2D}
end

(psum::ProfileSum)(x::Number, y::Number) = sum(prof -> prof(x, y), psum.profs)
intensity(psum::ProfileSum) = sum(intensity, psum.profs)
coord_mean(psum::ProfileSum) =
    sum(prof -> intensity(prof) * coord_mean(prof), psum.profs) / sum(intensity, psum.profs)

Base.:+(left::Profile2D, right::Profile2D) = ProfileSum{2}((left, right))
Base.:+(left::ProfileSum{N}, right::Profile2D) where {N} =
    ProfileSum{N + 1}((left.profs..., right))
Base.:+(left::Profile2D, right::ProfileSum{N}) where {N} =
    ProfileSum{N + 1}((left, right.profs...))
Base.:+(left::ProfileSum{N}, right::ProfileSum{M}) where {N,M} =
    ProfileSum{N + M}((left.profs..., right.profs...))

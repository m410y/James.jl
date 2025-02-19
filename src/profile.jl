abstract type Profile{N} end

const Profile2D = Profile{2}

function intensity(prof::Profile)
    solve(
        IntegralProblem((coord, p) -> prof(coord...), [-Inf, -Inf], [Inf, Inf]),
        HCubatureJL(),
    ).u
end

function coord_mean(prof::Profile)
    int =
        solve(
            IntegralProblem(
                (coord, p) -> [1, coord...] * prof(coord...),
                [-Inf, -Inf],
                [Inf, Inf],
            ),
            HCubatureJL(),
        ).u
    int[2:end] / int[1]
end

function coord_disp(prof::Profile)
    int =
        solve(
            IntegralProblem(
                (coord, p) -> [1, coord...] * [1, coord...]' * prof(coord...),
                [-Inf, -Inf],
                [Inf, Inf],
            ),
            HCubatureJL(),
        ).u
    (int[2:end, 2:end] * int[1, 1] - int[2:end, 1] * int[1, 2:end]') / int[1, 1]^2
end

struct Gaussian{N} <: Profile{N} end
(::Gaussian{N})(coord::Vararg{Number,N}) where {N} =
    exp(-sum(abs2, coord) / 2) * (2π)^(-N / 2)
intensity(::Gaussian) = 1
coord_mean(::Gaussian{N}) where {N} = Vec(zeros(N))
coord_disp(::Gaussian) = Diagonal(ones(N))

struct GriddedProfile{N} <: Profile{N}
    ext::Interpolations.AbstractExtrapolation
end

(prof::GriddedProfile{N})(coord::Vararg{Number,N}) where {N} = prof.ext(coord...)

struct AffinedProfile{N,P<:Profile{N}} <: Profile{N}
    profile::P
    A::Number
    V::Vec{N}
    M::Mat{N,N}
end

AffinedProfile{N}(prof::Profile{N}, args::Vararg{Number}) where {N} =
    AffinedProfile(p[1], args[1], Vec{N}(args[2:1+N]), Mat{Tuple{N,N}}(args[2+N:end]))

function (prof::AffinedProfile{N})(coord::Vararg{Number,N}) where {N}
    new_coord = prof.M * collect(coord) + prof.V
    prof.A * prof.profile(new_coord...)
end
intensity(prof::AffinedProfile) = prof.A / det(prof.M) * intensity(prof.profile)
coord_mean(prof::AffinedProfile) = inv(prof.M) * (coord_mean(prof.profile) - prof.V)
coord_disp(prof::AffinedProfile) = inv(prof.M) * coord_disp(prof.profile) * inv(prof.M)'

struct ProfileSum{N,M} <: Profile{N}
    profs::NTuple{M,Profile{N}}
end

(prof::ProfileSum)(coord::Vararg{Number,N}) where {N} = sum(p -> p(coord...), prof.profs)
intensity(prof::ProfileSum) = sum(intensity, prof.profs)
coord_mean(prof::ProfileSum) =
    sum(p -> intensity(p) * coord_mean(p), prof.profs) / sum(intensity, prof.profs)

Base.:+(left::Profile{N}, right::Profile{N}) where {N} = ProfileSum{N,2}((left, right))
Base.:+(left::ProfileSum{N,M}, right::Profile{N}) where {N,M} =
    ProfileSum{N,M + 1}((left.profs..., right))
Base.:+(left::Profile{N}, right::ProfileSum{N,M}) where {N,M} =
    ProfileSum{N,M + 1}((left, right.profs...))
Base.:+(left::ProfileSum{N,M1}, right::ProfileSum{N,M2}) where {N,M1,M2} =
    ProfileSum{N,M1 + M2}((left.profs..., right.profs...))

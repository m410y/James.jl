abstract type Goniometer{N} end

struct MultiaxisGoniometer{N} <: Goniometer{N}
    axes::NTuple{N,Axis}
    prelim::Transformation
end

const OneCircleGoniometer = MultiaxisGoniometer{1}
const TwoCircleGoniometer = MultiaxisGoniometer{2}
const ThreeCircleGoniometer = MultiaxisGoniometer{3}

function MultiaxisGoniometer(axes::Vararg{Union{Transformation,Axis}}; prelim = IdentityTransformation())
    active_axes = Axis[]
    for axis in axes[end:-1:begin]
        if axis isa Transformation
            prelim = prelim ∘ axis
        else
            push!(active_axes, prelim(axis))
        end
    end
    MultiaxisGoniometer(Tuple(active_axes[end:-1:begin]), prelim)
end

function fix_angle(gonio::MultiaxisGoniometer, (n, angle)::Pair)
    axes = Any[gonio.axes...]
    axes[n] = axes[n](angle)
    MultiaxisGoniometer(axes..., prelim = gonio.prelim)
end

function (gonio::MultiaxisGoniometer{N})(angles::Vararg{Number,N}) where {N}
    trans = gonio.prelim
    for (angle, axis) in zip(angles, gonio.axes)
        trans = axis(angle) ∘ trans
    end
    trans
end

function orient_angles(gonio::TwoCircleGoniometer, src::AbstractVector, dst::AbstractVector)
    axis₁, axis₂ = (axis.v for axis in gonio.axes)
    S = normalize(gonio.prelim(src))
    D = normalize(dst)
    orient_angles(axis₁, axis₂, S, D)
end

function reflection_angles(gonio::OneCircleGoniometer, s::AbstractVector, k::AbstractVector)
    A = gonio.axes[1].v
    S = normalize(gonio.prelim(s))
    K = normalize(k)
    θ = asin(norm(s) / 2norm(k))u"rad"
    reflection_angles(A, S, K, θ)
end

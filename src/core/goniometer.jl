struct Goniometer{N}
    axes::NTuple{N,Axis}
    prelim::Transformation
end

const OneCircleGoniometer = Goniometer{1}
const TwoCircleGoniometer = Goniometer{2}
const ThreeCircleGoniometer = Goniometer{3}

function Goniometer(
    axes::Vararg{Union{Transformation,Axis}};
    prelim = IdentityTransformation(),
)
    active_axes = Axis[]
    for axis in axes[end:-1:begin]
        if axis isa Transformation
            prelim = prelim ∘ axis
        else
            push!(active_axes, prelim(axis))
        end
    end
    Goniometer(Tuple(active_axes[end:-1:begin]), prelim)
end

function fix_angle(gonio::Goniometer, (n, angle)::Pair)
    axes = Any[gonio.axes...]
    axes[n] = axes[n](angle)
    Goniometer(axes..., prelim = gonio.prelim)
end

function (gonio::Goniometer{N})(angles::Vararg{Number,N}) where {N}
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

function Base.show(io::IO, ::MIME"text/plain", gonio::Goniometer{N}) where {N}
    print(io, "$N-circle Goniometer:\n")
    for (n, axis) in enumerate(gonio.axes)
        print(io, "  axis $n:\n")
        print(io, "    direction: [$(axis.v[1]), $(axis.v[2]), $(axis.v[3])]\n")
        print(io, "    position: [$(axis.p[1]), $(axis.p[2]), $(axis.p[3])]\n")
    end
    print(io, "  preliminary transform:\n")
    print(io, "    $(gonio.prelim)")
end

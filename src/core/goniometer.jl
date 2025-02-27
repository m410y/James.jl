struct Goniometer{N,T}
    axes::NTuple{N,Axis{T}}
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
    for axis in reverse(axes)
        if axis isa Transformation
            prelim = prelim ∘ axis
        else
            push!(active_axes, prelim(axis))
        end
    end
    Goniometer(Tuple(reverse(active_axes)), prelim)
end

function fix_angle(gonio::Goniometer, (n, angle))
    axes = Any[gonio.axes...]
    axes[n] = axes[n](angle)
    Goniometer(axes..., prelim = gonio.prelim)
end

function (gonio::Goniometer)(angles...)
    trans = gonio.prelim
    for (angle, axis) in zip(angles, gonio.axes)
        trans = axis(angle) ∘ trans
    end
    trans
end

function Base.show(io::IO, ::MIME"text/plain", gonio::Goniometer{N}) where {N}
    println(io, "$N-circle Goniometer:")
    for (n, axis) in enumerate(gonio.axes)
        p = axis.p * SpaceUnit |> u"μm"
        println(io, "  axis $n:")
        println(io, @sprintf("    position: [%6.2f%6.2f%6.2f]", p...))
        println(io, @sprintf("    direction: [%6.3f%6.3f%6.3f]", axis.v...))
    end
end

struct Motorized{G<:Goniometer,O}
    goniometer::G
    object::O
end

function (motor::Motorized)(angles...)
    motor.goniometer(angles...)(motor.object)
end

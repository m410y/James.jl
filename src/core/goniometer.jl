struct Goniometer{N,T<:Transformation}
    axes::NTuple{N,Axis}
    prelim::T
end

const OneCircleGoniometer = Goniometer{1}
const TwoCircleGoniometer = Goniometer{2}
const ThreeCircleGoniometer = Goniometer{3}

# TODO: optimize
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

# TODO: optimize
function fix_angle(gonio::Goniometer, (n, angle))
    axes = Any[gonio.axes...]
    axes[n] = axes[n](angle)
    Goniometer(axes..., prelim = gonio.prelim)
end

function (gonio::Goniometer{N})(angles::Vararg{Number,N}) where {N}
    trans = gonio.prelim
    for n = 1:N
        trans = gonio.axes[n](angles[n]) ∘ trans
    end
    trans
end

function Base.show(io::IO, ::MIME"text/plain", gonio::Goniometer{N}) where {N}
    punit = u"μm"
    println(io, summary(gonio), ":")
    for (n, axis) in enumerate(gonio.axes)
        println(io, "  Axis $n:")
        println(
            io,
            "    position [$punit]: ",
            @sprintf("%6.2f, %6.2f, %6.2f", space_convert.(punit, axis.p)...)
        )
        println(io, "    direction    : ", @sprintf("%6.3f, %6.3f, %6.3f", axis.v...))
    end
end

struct Motorized{G<:Goniometer,O}
    goniometer::G
    object::O
end

function (motor::Motorized)(angles...)
    motor.goniometer(angles...)(motor.object)
end

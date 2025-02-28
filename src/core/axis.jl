struct Axis{T}
    v::Vec3{T}
    p::Point3{T}
    function Axis{T}(v, p) where {T}
        n = normalize(v)
        pn = p - n * dot(n, p)
        new{T}(n, pn)
    end
end

Axis(v, p) = Axis{Float64}(v, p)

function (axis::Axis)(angle)
    rot = AngleAxis(angle, axis.v..., false)
    AffineMap(rot, axis.p - rot * axis.p)
end

(::IdentityTransformation)(axis::Axis) = axis
function (trans::LinearMap)(axis::Axis)
    Axis(trans(axis.v), trans(axis.p))
end
function (trans::Translation)(axis::Axis)
    Axis(axis.v, trans(axis.p))
end
function (trans::AffineMap)(axis::Axis)
    Axis(trans.linear * axis.v, trans(axis.p))
end

function Base.show(io::IO, ::MIME"text/plain", axis::Axis; unit = u"μm")
    println(io, summary(axis), ":")
    println(io, "  position [$unit]: ", @sprintf("%6.2f, %6.2f, %6.2f", NoUnits.(axis.p * SpaceUnit / unit)...))
    println(io, "  direction    : ", @sprintf("%6.3f, %6.3f, %6.3f", axis.v...))
end

function orient_angles(axis₁::Axis, axis₂::Axis, src, dst)
    S = normalize(src)
    D = normalize(dst)
    orient_angles(axis₁.v, axis₂.v, S, D)
end

function reflection_angles(axis::Axis, s, k)
    S = normalize(s)
    K = normalize(k)
    θ = asin(norm(s) / 2norm(k))
    reflection_angles(axis.v, S, K, θ)
end

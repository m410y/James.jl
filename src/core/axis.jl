struct Axis
    v::Vec3
    p::Point3
    Axis(v, p) = new(normalize(v), p)
end

function (axis::Axis)(angle::Number)
    rot = AngleAxis(angle, axis.v..., false)
    AffineMap(rot, axis.p - rot * axis.p)
end

(::IdentityTransformation)(axis::Axis) = axis
function (trans::LinearMap)(axis::Axis)
    Axis(trans(axis.v), trans(axis.p))
end
function (trans::Translation)(axis::Axis)
    Axis(trans(axis.v), trans(axis.p))
end
function (trans::AffineMap)(axis::Axis)
    Axis(trans(axis.v), trans(axis.p))
end

function Base.show(io::IO, ::MIME"text/plain", axis::Axis)
    print(io, "Axis:\n")
    print(io, "  direction: [$(axis.v[1]), $(axis.v[2]), $(axis.v[3])]\n")
    print(io, "  position: [$(axis.p[1]), $(axis.p[2]), $(axis.p[3])]")
end

struct Axis
    p::Point3
    v::Vec3
    Axis(p, v) = new(p, normalize(v))
end

function (axis::Axis)(angle::Number)
    rot = AngleAxis(angle, axis.v..., false)
    AffineMap(rot, axis.p - rot * axis.p)
end

(::IdentityTransformation)(axis::Axis) = copy(axis)
function (trans::LinearMap)(axis::Axis)
    Axis(trans(axis.p), trans(axis.v))
end
function (trans::Translation)(axis::Axis)
    Axis(trans(axis.p), trans(axis.v))
end
function (trans::AffineMap)(axis::Axis)
    Axis(trans(axis.p), trans(axis.v))
end

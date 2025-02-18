struct Detector2D
    size::NTuple{2,Number}
    p::Point3
    ex::Vec3
    ey::Vec3
    axis::Axis
end

function rotate(detector::Detector2D, angle::Number)
    rot = detector.axis(angle)
    Detector2D(
        detector.size,
        rot(detector.p),
        rot(detector.ex),
        rot(detector.ey),
        detector.axis,
    )
end

(::IdentityTransformation)(detector::Detector2D) = copy(detector)
function (trans::LinearMap)(detector::Detector2D)
    Detector2D(
        detector.size,
        trans(detector.p),
        trans(detector.ex),
        trans(detector.ey),
        detector.axis,
    )
end
function (trans::Translation)(detector::Detector2D)
    Detector2D(
        detector.size,
        trans(detector.p),
        trans(detector.ex),
        trans(detector.ey),
        detector.axis,
    )
end
function (trans::AffineMap)(detector::Detector2D)
    Detector2D(
        detector.size,
        trans(detector.p),
        trans(detector.ex),
        trans(detector.ey),
        detector.axis,
    )
end

function (detector::Detector2D)(x::Number, y::Number)
    detector.p + detector.ex * x + detector.ey * y
end

intersect(detector::Detector2D, p::AbstractVector, v::AbstractVector) =
    (p - detector.p) \ [ex ey -v]

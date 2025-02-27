struct Detector{N,T}
    size::NTuple{N,Int}
    p::Point3{T}
    e::Mat{3,N,T}
end

const Detector2D = Detector{2}

(::IdentityTransformation)(detector::Detector) = detector
function (trans::LinearMap)(detector::Detector)
    Detector(detector.size, trans(detector.p), trans(detector.e))
end
function (trans::Translation)(detector::Detector)
    Detector(detector.size, trans(detector.p), detector.e)
end
function (trans::AffineMap)(detector::Detector)
    Detector(detector.size, trans(detector.p), trans.linear * detector.e)
end

(detector::Detector)(coord) = detector.p + detector.e * Vec(coord...)

function intersect_coord(detector::Detector, p, v)
    _, coord... = [-normalize(v) detector.e] \ (Vec(p...) - detector.p)
    Vec(coord...)
end

function Base.show(io::IO, ::MIME"text/plain", detector::Detector{N}) where {N}
    p = axis.p * SpaceUnit |> u"mm"
    println(io, "$N-dimentional Detector:")
    println(io, "  size: ", join(detector.size, "×"))
    println(io, @sprintf("  zero position: [%6.1f%6.1f%6.1f]", p...))
    println(io, "  coordinate lines:")
    for (n, col) in enumerate(eachcol(detector.e))
        e = col * SpaceUnit |> u"μm"
        println(io, @sprintf("    line %d: [%6.1f%6.1f%6.1f]", n, e...))
    end
end

struct Detector{N}
    size::NTuple{N,Number}
    p::Point3{Length}
    e::NTuple{N,Vec3{Length}}
end

const Detector2D = Detector{2}

(::IdentityTransformation)(detector::Detector) = detector
function (trans::LinearMap)(detector::Detector)
    Detector(detector.size, trans(detector.p), trans.(detector.e))
end
function (trans::Translation)(detector::Detector)
    Detector(detector.size, trans(detector.p), trans.(detector.e))
end
function (trans::AffineMap)(detector::Detector)
    Detector(detector.size, trans(detector.p), trans.(detector.e))
end

function (detector::Detector)(coord::Vararg{Number})
    detector.p + sum(detector.e .* coord)
end

function intersect_coord(detector::Detector, p::AbstractVector, v::AbstractVector)
    _, coord... = hcat(-ustrip(v), detector.e...) \ (p - detector.p)
    NoUnits.(coord)
end

function Base.show(io::IO, ::MIME"text/plain", detector::Detector)
    print(io, "$(length(detector.size))-dimentional Detector:\n")
    print(io, "  size: $(detector.size)\n")
    print(io, "  zero position: [$(detector.p[1]), $(detector.p[2]), $(detector.p[3])]\n")
    print(io, "  coordinate lines:\n")
    for (n, e) in enumerate(detector.e)
        print(io, "    line $n: [$(e[1]), $(e[2]), $(e[3])]\n")
    end
end

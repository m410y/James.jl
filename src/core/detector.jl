struct Detector{N,T<:AbstractAffineMap}
    size::NTuple{N,Int}
    trans::T
end

const Detector2D = Detector{2}

(::IdentityTransformation)(detector::Detector) = detector
function (trans::LinearMap)(detector::Detector)
    Detector(detector.size, trans ∘ detector.trans)
end
function (trans::Translation)(detector::Detector)
    Detector(detector.size, trans ∘ detector.trans)
end
function (trans::AffineMap)(detector::Detector)
    Detector(detector.size, trans ∘ detector.trans)
end

(detector::Detector{N})(coord::Vararg{Number,N}) where {N} = detector.trans(Vec{N}(coord))
(detector::Detector{N})(coord::AbstractVector) where {N} = detector.trans(coord)

# TODO: make it work for arbitrary transform
function intersect_coord(detector::Detector{2}, p, v)
    _, coord... =
        [-Vec3(v) detector.trans.linear] \ (Vec3(p...) - detector.trans.translation)
    Vec2(coord...)
end

function Base.show(io::IO, ::MIME"text/plain", detector::Detector{N}) where {N}
    punit = u"mm"
    eunit = u"μm"
    p = detector(zeros(N))
    println(io, summary(detector), ":")
    println(io, "  size: ", join(detector.size, "×"))
    println(
        io,
        "  zero position [$punit]: ",
        @sprintf("%6.1f, %6.1f, %6.1f", space_convert.(punit, p)...)
    )
    println(io, "  coordinate lines [$eunit]:")
    for n = 1:N
        coord = zeros(N)
        coord[n] = 1
        e = detector(coord) - p
        println(
            io,
            "    line $n: ",
            @sprintf("%6.1f, %6.1f, %6.1f", space_convert.(eunit, e)...)
        )
    end
end

abstract type Sample end

struct SingleCrystal{T} <: Sample
    p::Point3{T}
    UB::Mat3{T}
end

(trans::IdentityTransformation)(cryst::SingleCrystal) = cryst
function (trans::LinearMap)(cryst::SingleCrystal)
    SingleCrystal(trans(cryst.p), trans(cryst.UB))
end
function (trans::Translation)(cryst::SingleCrystal)
    SingleCrystal(trans(cryst.p), cryst.UB)
end
function (trans::AffineMap)(cryst::SingleCrystal)
    SingleCrystal(trans(cryst.p), trans.linear * cryst.UB)
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    cryst::SingleCrystal;
    punit = u"μm",
    runit = u"eV",
)
    println(io, summary(cryst), ":")
    println(
        io,
        "  position [$punit]: ",
        @sprintf("%6.1f, %6.1f, %6.1f", NoUnits.(cryst.p * SpaceUnit / punit)...)
    )
    println(io, "  orientation matrix (UB) [$runit]:")
    convert(E) =
        runit != ReciprocalUnit ? ustrip(uconvert(runit, E * ReciprocalUnit, Spectral())) :
        E
    for row in eachrow(convert.(cryst.UB))
        println(io, @sprintf("    %6.0f %6.0f %6.0f", row...))
    end
end

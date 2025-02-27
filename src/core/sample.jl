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

function Base.show(io::IO, ::MIME"text/plain", cryst::SingleCrystal)
    p = axis.p * SpaceUnit |> u"μm"
    UB = round.(cryst.UB, sigdigits = sigdigits)
    println(io, "SingleCrystal:")
    println(io, @sprintf("  position: [%6.1f%6.1f%6.1f]", p...))
    println(io, "  orientation matrix (UB) in $(ReciprocalUnit):")
    for row in eachrow(UB)
        println(io, @sprintf("    %6.0f%6.0f%6.0f", row...))
    end
end

abstract type Crystal end

struct SingleCrystal <: Crystal
    pos::Point3
    UB::Mat3
end

(trans::IdentityTransformation)(cryst::SingleCrystal) = cryst
function (trans::LinearMap)(cryst::SingleCrystal)
    SingleCrystal(trans(cryst.pos), trans.linear * cryst.UB)
end
function (trans::Translation)(cryst::SingleCrystal)
    SingleCrystal(trans(cryst.pos), trans.linear * cryst.UB)
end
function (trans::AffineMap)(cryst::SingleCrystal)
    SingleCrystal(trans(cryst.pos), trans.linear * cryst.UB)
end

function Base.show(io::IO, ::MIME"text/plain", cryst::SingleCrystal)
    UBs = round.(ustrip.(cryst.UB), sigdigits = 5)
    print(io, "SingleCrystal:\n")
    print(io, "  position: [$(cryst.pos[1]), $(cryst.pos[2]), $(cryst.pos[3])]\n")
    print(io, "  orientation matrix (UB) in $(Unitful.unit(first(cryst.UB))):\n")
    print(io, "    $(UBs[1, 1]) $(UBs[1, 2]) $(UBs[1, 3])\n")
    print(io, "    $(UBs[2, 1]) $(UBs[2, 2]) $(UBs[2, 3])\n")
    print(io, "    $(UBs[3, 1]) $(UBs[3, 2]) $(UBs[3, 3])")
end

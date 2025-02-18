abstract type Crystal end

struct SingleCrystal <: Crystal
    pos::Point3
    UB::Mat3
end

reflex(cryst::SingleCrystal, hkl::AbstractVector) = Vec3(cryst.UB * hkl)

(trans::IdentityTransformation)(cryst::SingleCrystal) = copy(cryst)
function (trans::LinearMap)(cryst::SingleCrystal)
    SingleCrystal(trans(cryst.pos), trans.linear * cryst.UB)
end
function (trans::Translation)(cryst::SingleCrystal)
    SingleCrystal(trans(cryst.pos), trans.linear * cryst.UB)
end
function (trans::AffineMap)(cryst::SingleCrystal)
    SingleCrystal(trans(cryst.pos), trans.linear * cryst.UB)
end

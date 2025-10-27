module D8Venture

using GeometryBasics
using Unitful
import Unitful: Length

using JamesCore
using James
using P4P

export collision

function angle_dif(d::Length)
    u"rad"(
        d < 159.357u"mm" ? 1.69095atan(61.47048u"mm"/d) :
        d < 209.205u"mm" ? evalpoly(d, (6.32162, -5.35416e-2u"mm^-1", 1.1149e-4u"mm^-2")) :
        0.0
    )
end

function collision(tthD, omega, distance)
   cos(tthD - omega + π/2) < cos(angle_dif(distance))
end

function angle_calculator(fname::AbstractString)
    axes = [RotAxis("-z"), RotAxis("-x"), RotAxis("z")]
    χ = 54.7112u"°"
    axes_reduced, prelim = fix_axes_params(axes, [2 => χ])
    n = Vec3("x")
    p4p = P4P.load(fname)
    UB = p4p["ORT"] 
    λ = p4p["SOURCE"][2]
    AngleCalculator(
        prelim.linear * UB,
        (ax.v for ax in axes_reduced),
        n / λ
    )
end

end

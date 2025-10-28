module D8Venture

using GeometryBasics
using CoordinateTransformations
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
    p4p = P4P.load(fname)
    UB = p4p["ORT"] 
    model = instrument_model(
        λ = p4p["SOURCE"][2]u"Å",
        xc = p4p["ADPAR"][1],
        yc = p4p["ADPAR"][2],
    )
    AngleCalculator(
        model.prelim * UB,
        (ax.v for ax in model.saxes),
        model.beam.k
    )
end

struct Model{T<:Real}
    daxes::Tuple{TransAxis{T},RotAxis{T}}
    saxes::Tuple{RotAxis{T},RotAxis{T}}
    prelim::Mat3{T}
    detect::Detector{T}
    beam::Ray{T}
end

function instrument_model(; 
    λ::Length = 1.54184u"Å",
    xc = 387.2198,
    yc = 503.9324,
    χ = 54.7112u"°",
    px::Length = 135.3u"μm"
)
    euler_axes = [RotAxis("-z"), RotAxis("-x"), RotAxis("z")]
    saxes, prelim = fix_axes_params(euler_axes, [2 => χ])
    prelim = Mat3(prelim.linear)
    k = Vec3("x") / NoUnits(λ / u"Å")
    ex = NoUnits(px / u"mm") * Vec3("-y")
    ey = NoUnits(px / u"mm") * Vec3("z")
    p = -xc * ex - yc * ey
    Model{Float64}(
        (TransAxis("x"), RotAxis("z")),
        tuple(saxes...),
        prelim,
        Detector([ex ey], p),
        Ray(k)
    )
end

end

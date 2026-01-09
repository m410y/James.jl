module D8Venture

using JamesCore: RigidObject
using GeometryBasics
using CoordinateTransformations
using Unitful
using LinearAlgebra
import Unitful: Length

using JamesCore
using James
using P4P
using SFRM

export collision, Model, angle_calculator, reflex_xy, fit_equator_reflex

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

struct Model{T<:Real}
    cryst::Lattice{T}
    daxes::Tuple{TransAxis{T},RotAxis{T}}
    saxes::NTuple{3,RotAxis{T}}
    detect::Plane{T}
    beams::NTuple{2,Ray{T}}
end

function Model{T}(filename::AbstractString; 
    px::Length = 135.3u"μm",
) where {T<:Real}
    p4p = P4P.load(filename)
    exy = NoUnits(px / u"mm") * [Vec3("-y") Vec3("z")]
    Model{T}(
        Lattice{T}(p4p["ORT"]),
        (TransAxis{T}("x"), RotAxis{T}("z")),
        (RotAxis{T}("-z"), RotAxis{T}("-x"), RotAxis{T}("z")),
        Plane{T}(exy, - exy * Vec2(p4p["ADPAR"][1:2]...)),
        (Ray{T}(Vec3("x") / p4p["SOURCE"][2]), Ray{T}(Vec3("x") / p4p["SOURCE"][3]))
    )
end

function angle_calculator(model::Model{T}; χ = 54.7112u"°") where {T<:Real}
    prelim, saxes = fix_axes(model.saxes, 2 => χ)
    AngleCalculator{T}(prelim.linear * model.cryst.UB, (ax.v for ax in saxes), model.beams[1].k)
end

function reflex_xy(model::Model{<:Real}, fname::AbstractString, hkl::AbstractVector)
    sfrm = SFRM.load(fname)
    strans = fix_axes(model.saxes, (sfrm.phi, sfrm.chi, sfrm.omega + sfrm.increment / 2))
    dtrans = fix_axes(model.daxes, (sfrm.distance * 10, sfrm.tth))
    ray = Ray(model.beams[1].k + strans(model.cryst.UB * hkl))
    detect = dtrans(model.detect)
    ray2xy(ray, detect)
end

function fit_equator_reflex(model::Model{<:Real}, fname::AbstractString, hkl::AbstractVector; kwargs...)
    sfrm = SFRM.load(fname)
    tthD = deg2rad(sfrm.tth)
    detect = fix_axes(model.daxes, (sfrm.distance * 10, tthD))(model.detect)
    d_inv = norm(model.cryst.UB * hkl)
    xys = map(model.beams) do beam
        ray = model.daxes[2](copysign(2asin(d_inv/2norm(beam.k)), sin(tthD)))(beam)
        ray2xy(ray, detect)
    end
    fit_peak(reverse(sfrm.image, dims = 1), xys...; kwargs...)
end

end

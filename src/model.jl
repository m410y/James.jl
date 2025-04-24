export Model

const wavelengths = Dict(
    "Cu" => (Kα1 = 1.54059290u"Å", Kα2 = 1.54442740u"Å"),
    "Mo" => (Kα1 = 0.70931715u"Å", Kα2 = 0.71360700u"Å")
)

struct Crystal
    a::PVec
    b::PVec
    c::PVec
    p::PVec
end

struct Model
    beam::PVec
    E1::Real
    E2::Real
    detector::Detector
    phi_ax::BiVec
    omega_ax::BiVec
    chi_ax::BiVec
    tth_ax::BiVec
end

function model_from_sfrm(sfrm::SiemensFrame; correct_wavelen = false, px = missing, center = missing, distance = missing)
    λ1, λ2 = if correct_wavelen
        anode = titlecase(sfrm.target)
        wavelengths[anode].Kα1, wavelengths[anode].Kα2
    else
        sfrm.lambdaKα1 * u"Å", sfrm.lambdaKα2 * u"Å"
    end
    E1 = ustrip(uconvert(runit, λ1, Spectral()))
    E2 = ustrip(uconvert(runit, λ2, Spectral()))
    distance = coalesce(distance, ustrip(uconvert(munit, sfrm.distance * u"cm")))
    px = coalesce(px, ustrip(uconvert(munit, 512u"cm" / sfrm.pix512percm / size(sfrm.image, 1))))
    cx, cy = coalesce(center, (sfrm.xcenter, sfrm.ycenter))
    ex = direction(0, -px, 0)
    ey = direction(0, 0, px)
    pcenter = point(distance, 0, 0)
    Model(
        direction(1, 0 ,0),
        E1,
        E2,
        detector = Detector(
            ex,
            ey,
            pcenter - cx * ex - cy * ey
        ),
        rotaxis(0, 0, -1),
        rotaxis(0, 0, 1),
        rotaxis(-1, 0, 0),
        rotaxis(0, 0, 1),
    )
end
model_from_sfrm(path::AbstractString; kwargs...) = model_from_sfrm(SFRM.load(path), kwargs...)

function crystal_from_p4p(p4p::AbstractDict; pos = missing)
    UB = p4p["ORT"] * ustrip(uconvert(runit, 1u"Å^-1", Spectral()))
    p = coalesce(pos, zeros(3))
    Crystal(
        (direction(col) for col in eachcol(UB))...,
        point(p...)
    )
end
crystal_from_p4p(path::AbstractString; kwargs...) = crystal_from_p4p(P4P.load(path), kwargs...)

struct Frame{T<:Real}
    image::AbstractMatrix{T}
    tth::Real,
    phi::Real
    chi::Real
    omega::Real
    axis::Symbol
    inc::Real
end

function frame_from_sfrm(sfrm::SiemensFrame)
    Frame(
        sfrm.image,
        sfrm.tth,
        sfrm.phi,
        sfrm.chi,
        sfrm.omega,
        (:tth, :omega, :phi, :chi)[sfrm.axis],
        sfrm.increment
    )
end

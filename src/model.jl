export Model

struct Model
    xray1::XRay
    xray2::XRay
    detector::Detector
    theta_axis::RotAxis
    omega_axis::RotAxis
    chi_axis::RotAxis
    phi_axis::RotAxis
end

function Model(sfrm::SiemensFrame)
    E1 = uconvert(runit, sfrm.lambdaKα1 * u"Å", Spectral()) |> ustrip
    E2 = uconvert(runit, sfrm.lambdaKα2 * u"Å", Spectral()) |> ustrip
    distance = uconvert(munit, sfrm.distance * u"cm") |> ustrip
    px = uconvert(munit, 512u"cm" / sfrm.pix512percm) |> ustrip
    center = [sfrm.xcenter, sfrm.ycenter]
    Model(
        XRay(E1, 0, 0),
        XRay(E2, 0, 0),
        Detector([0, -px, 0], [0, 0, px], [distance, center[1] * px, -center[2] * px]),
        RotAxis(0, 0, 1),
        RotAxis(0, 0, 1),
        RotAxis(-1, 0, 0),
        RotAxis(0, 0, -1)
    )
end
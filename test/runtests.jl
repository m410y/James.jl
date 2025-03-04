using Test
using James
using FileIO
using LinearAlgebra
using Rotations
using GeometryBasics
using CoordinateTransformations
using Unitful, UnitfulEquivalences
using Distributions

# include("Aqua.jl")

function Kα_spec(amplitude::Real, λ::Unitful.Length, FWHM::Unitful.Energy, div::Real)
    E = uconvert(u"eV", λ, Spectral()) |> ustrip
    ΔE = uconvert(u"eV", FWHM) / 2 |> ustrip
    Δk = E * div
    dist = product_distribution(
        Cauchy(E, ΔE),
        Distributions.Normal(0, Δk),
        Distributions.Normal(0, Δk),
    )
    mean = [E, 0, 0]
    intensity = amplitude / pdf(dist, mean)
    cov = [
        4ΔE^2 0 0
        0 Δk^2 0
        0 0 Δk^2
    ]
    James.DistSpectrum(intensity, mean, cov, dist)
end

@testset "James" begin
    ϕ_axis = James.Axis([0, 0, -1], zero(Point3))
    @test ϕ_axis(1.2345).linear ≈ Matrix(RotZ(-1.2345)) rtol = 1e-15
    χ_axis = James.Axis([-1, 0, 0], zero(Point3))
    @test χ_axis(1.2345).linear ≈ Matrix(RotX(-1.2345)) rtol = 1e-15
    ω_axis = James.Axis([0, 0, 1], zero(Point3))
    @test ω_axis(1.2345).linear ≈ Matrix(RotZ(1.2345)) rtol = 1e-15
    gonio_euler = James.Goniometer(ϕ_axis, χ_axis, ω_axis)
    gonio_sample = James.fix_angle(gonio_euler, 2 => deg2rad(54.7112))
    gonio_detector = James.Goniometer(James.Axis([0, 0, 1], zero(Point3)))
    detector_0 = let center_coord = [388.88, 504.83], pixel_size = 0.1353
        linear = Mat{3,2}([
            0 0
            -pixel_size 0
            0 pixel_size
        ])
        translation = Vec3(-linear * center_coord)
        James.Detector((768, 1024), AffineMap(linear, translation))
    end
    @test detector_0(388.88, 504.83) ≈ zero(Point3) rtol = 1e-15
    MoKα1 = Kα_spec(2e4, 0.70931715u"Å", 6.31u"eV", 1e-3)
    MoKα2 = Kα_spec(1e4, 0.713607u"Å", 6.49u"eV", 1e-3)
    MoKα = James.SpectrumSum(MoKα1, MoKα2)
    @test mean(MoKα) ≈ [17450, 0, 0] rtol = 1e-3
    CuKα1 = Kα_spec(2e4, 1.54059290u"Å", 2.11u"eV", 1e-3)
    CuKα2 = Kα_spec(1e4, 1.54442740u"Å", 2.17u"eV", 1e-3)
    CuKα = James.SpectrumSum(CuKα1, CuKα2)
    @test mean(CuKα) ≈ [8048, 0, 0] rtol = 1e-3
    preset = Dict(
        "spectrum_cu" => CuKα,
        "spectrum_mo" => MoKα,
        "detector" => detector_0,
        "detector_goniometer" => gonio_detector,
        "sample_goniometer" => gonio_sample,
    )
    p4p = load("data/20240610_Ge_on_XiHead_5.p4p")
    sfrm = load(
        "data/20240703_Ge_on_Chi_bond_2/mo_Ge_on_Chi_5_10_0_6_family_M93p9_phiM178p60_01_0001.sfrm",
    )
    crystal = James.SingleCrystal(
        zero(Point3),
        Mat3([
            263.8191691135128 -29.23875709212313 -2173.4789539515464
            -1251.6637536684098 1787.9694304584657 -175.9809519389325
            1777.1334594838004 1263.6362877901117 198.71123540450782
        ]),
    )
    crystal_loaded = James.load_crystal_p4p(p4p)
    @test crystal_loaded.UB ≈ crystal.UB rtol = 1e-10
    frame = James.load_frame_sfrm(sfrm)
    setting = James.ScanFrameSetting(
        Vec1(deg2rad(266.1)),
        Vec2(deg2rad(181.4004), deg2rad(245.62)),
        2,
        deg2rad(4)
    )
    @test frame.setting == setting
    experiment = James.load_experiment(sfrm, p4p, preset)
    @test experiment.detector.object.trans.translation[1] ≈ 128.533 atol = 1e-3
    save("data/experiment.jld2", "exp", experiment)
end

function load_experiment(sfrm, p4p, preset)
    spectrum = preset["spectrum_$(lowercase(sfrm["TARGET"]))"]
    detector = Translation([sfrm["DISTANC"][2], 0, 0]u"cm")(preset["detector"])
    gonio_d = preset["detector_goniometer"]
    crystal = SingleCrystal(zeros(3)u"m", uconvert.(u"eV", p4p["ORT"]u"Å^-1", Spectral()))
    gonio_c = preset["sample_goniometer"]
    Experiment(spectrum, Motorized(gonio_d, detector), Motorized(gonio_c, crystal))
end

function load_frame(sfrm)
    image = sfrm["IMG"]
    angles_d = (sfrm["ANGLES"][1]u"°",)
    angles_s = (sfrm["ANGLES"][3]u"°", sfrm["ANGLES"][2]u"°")
    setting = if sfrm["INCREME"] == 0
        StillFrameSetting(angles_d, angles_s)
    else
        axis =
            sfrm["AXIS"] == 2 ? 2 :
            sfrm["AXIS"] == 3 ? 1 : error("not supported axis $(sfrm["AXIS"])")
        increment = sfrm["INCREME"]u"°"
        ScanFrameSetting(Vec(angles_d), Vec(angles_s), axis, increment)
    end
    Frame(image, setting)
end

function load_crystal_p4p(p4p_file::AbstractDict; unit = u"Å^-1")
    coef = ReciprocalUnit != unit ? ustrip(uconvert(ReciprocalUnit, 1.0 * unit, Spectral())) : 1.0
    SingleCrystal(zero(Point3), Mat3(p4p_file["ORT"]) * coef)
end

load_crystal_p4p(filename::AbstractString; unit = u"Å^-1") = load_crystal_p4p(load(filename), unit = unit)

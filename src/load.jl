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
            sfrm["AXIS"] == 2 ? model.gonio_c.axes[2] :
            sfrm["AXIS"] == 3 ? model.gonio_c.axes[1] :
            error("not supported axis $(sfrm["AXIS"])")
        increment = sfrm["INCREME"]u"°"
        ScanFrameSetting(angles_d, angles_s, axis, increment)
    end
    Frame(image, setting)
end

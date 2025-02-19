function load_model(sfrm, p4p, preset)
    spectrum = preset["spectrum_$(lowercase(sfrm["TARGET"]))"]
    detector = let
        translation = Translation([sfrm["DISTANC"][2], 0, 0]u"cm")
        rotation = preset["detector_goniometer"](sfrm["ANGLES"][1]u"°")
        detector_0 = preset["detector"]
        rotation(translation(detector_0))
    end
    crystal = let
        rotation = preset["sample_goniometer"](sfrm["ANGLES"][3]u"°", sfrm["ANGLES"][2]u"°")
        crystal_0 =
            SingleCrystal(zeros(3)u"m", uconvert.(u"eV", p4p["ORT"]u"Å^-1", Spectral()))
        rotation(crystal_0)
    end
    if sfrm["AXIS"] == 2
        axis = preset["sample_goniometer"].axes[2]
        angles = (0u"°", sfrm["INCREME"]u"°")
    elseif sfrm["AXIS"] == 3
        axis = preset["sample_goniometer"].axes[1]
        angles = (0u"°", sfrm["INCREME"]u"°")
    else
        error("not supported axis $(sfrm["AXIS"])")
    end
    ScanModel(spectrum, detector, crystal, axis, angles)
end

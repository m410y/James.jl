function load_moving_model(sfrm::AbstractDict, p4p::AbstractDict, preset::AbstractDict)
    spectrum = preset["spectrum_$(lowercase(sfrm["TARGET"]))"]
    detector = Translation([sfrm["DISTANC"][2], 0, 0]u"cm")(preset["detector"])
    gonio_d = preset["detector_goniometer"]
    crystal = SingleCrystal(zeros(3)u"m", uconvert.(u"eV", p4p["ORT"]u"Å^-1", Spectral()))
    gonio_c = preset["sample_goniometer"]
    MovingModel(spectrum, detector, gonio_d, crystal, gonio_c)
end

function load_still_model(sfrm::AbstractDict, p4p::AbstractDict, preset::AbstractDict)
    model = load_moving_model(sfrm, p4p, preset)
    detector = model.gonio_d(sfrm["ANGLES"][1]u"°")(model.detect)
    crystal = model.gonio_c(sfrm["ANGLES"][3]u"°", sfrm["ANGLES"][2]u"°")(model.cryst)
    StillModel(model.spec, detector, crystal)
end

function load_scan_model(sfrm::AbstractDict, p4p::AbstractDict, preset::AbstractDict)
    model = load_still_model(sfrm, p4p, preset)
    axis =  sfrm["AXIS"] == 2 ? model.gonio_c.axes[2] :
            sfrm["AXIS"] == 3 ? model.gonio_c.axes[1] :
            error("not supported axis $(sfrm["AXIS"])")
    angles = (0u"°", sfrm["INCREME"]u"°")
    ScanModel(model.spec, model.detect, model.cryst, axis, angles)
end

load_image(sfrm::AbstractDict) = sfrm["IMG"]

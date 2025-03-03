function load_experiment(sfrm::AbstractDict, p4p::AbstractDict, preset::AbstractDict)
    spectrum = preset["spectrum_$(lowercase(sfrm["TARGET"]))"]
    distanc_coef = ReciprocalUnit != u"cm" ? ustrip(uconvert(SpaceUnit, 1u"cm")) : 1
    detector = Translation([sfrm["DISTANC"][2], 0, 0] * distanc_coef)(preset["detector"])
    gonio_d = preset["detector_goniometer"]
    crystal = load_crystal_p4p(p4p)
    gonio_c = preset["sample_goniometer"]
    Experiment(spectrum, Motorized(gonio_d, detector), Motorized(gonio_c, crystal))
end

function load_frame_sfrm(sfrm::AbstractDict)
    image = sfrm["IMG"]
    angles_d = (sfrm["ANGLES"][1],) .|> deg2rad
    angles_s = (sfrm["ANGLES"][3], sfrm["ANGLES"][2]) .|> deg2rad
    setting = if sfrm["INCREME"] == 0
        StillFrameSetting(angles_d, angles_s)
    else
        axis =
            sfrm["AXIS"] == 2 ? 2 :
            sfrm["AXIS"] == 3 ? 1 : error("not supported axis $(sfrm["AXIS"])")
        increment = sfrm["INCREME"] |> deg2rad
        ScanFrameSetting(Vec(angles_d), Vec(angles_s), axis, increment)
    end
    Frame(image, setting)
end

load_frame_sfrm(filename::AbstractString) = load_frame_sfrm(load(filename))

function load_crystal_p4p(p4p_file::AbstractDict)
    coef =
        ReciprocalUnit != u"Å^-1" ? ustrip(uconvert(ReciprocalUnit, 1u"Å^-1", Spectral())) :
        1
    SingleCrystal(zero(Point3), Mat3(p4p_file["ORT"]) * coef)
end

load_crystal_p4p(filename::AbstractString) = load_crystal_p4p(load(filename))

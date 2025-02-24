function describe_peak(
    img::AbstractArray,
    indices::AbstractArray{CartesianIndex{N}},
) where {N}
    coord = sum(idx -> SVector{N}(idx.I...) * img[idx], indices) / sum(img[indices])
    size = sqrt(length(indices))
    amplitude = maximum(img[indices])
    coord, size, amplitude
end

function collect_peaks_higher(img::AbstractArray; baseline::Number = 2median(img))
    labeled_mask = label_components(img .> baseline)
    map(
        i -> describe_peak(img, i),
        component_indices(CartesianIndex, labeled_mask)[begin+1:end],
    )
end

function reflex(model::StillModel, coord)
    r = model.detect(coord...) - model.cryst.pos
    n = normalize(r)
    k0 = wvec_mean(model.spec)
    s = norm(k0) * n - k0
    Vec3(inv(Matrix(model.cryst.UB)) * s)
end

function estimate_reflex(
    wvec::AbstractVector,
    detector_gonio::Goniometer,
    detector::Detector,
    sample_gonio::Goniometer,
    detector_angle,
    sample_angle,
    coord;
    cryst_pos::Point3 = Point3(0, 0, 0)u"m",
)
    p = detector_gonio(detector_angle...)(detector)(coord...)
    s = Vec3(wvec - norm(wvec) * normalize(p - cryst_pos))
    inv(sample_gonio(sample_angle...))(s)
end

function combine_reflexes_fast(reflexes, ΔE)
    comb = Dict()
    for s in reflexes
        rs = round.(s ./ ΔE)
        if !haskey(comb, rs)
            comb[rs] = [s]
        else
            push!(comb[rs], s)
        end
    end
    mean.(values(comb))
end

function combine_reflexes_precise(reflexes, max_dif)
    comb = [[zero(first(reflexes))]]
    for s in reflexes
        dif, idx = findmin(s_arr -> norm(first(s_arr) - s), comb)
        if dif > max_dif
            push!(comb, [s])
        else
            push!(comb[idx], s)
        end
    end
    popfirst!(comb)
    map(mean, comb)
end

function collect_sfrm_orient(files)
    df = DataFrame(:angles => [], :peak => [])
    for file in files
        if all(splitext(file)[2] .!= FileIO.info(format"SFRM")[2])
            continue
        end
        sfrm = try
            load(file)
        catch
            println("$(file): failed")
            continue
        end
        angles = (sfrm["ANGLES"] + sfrm["ENDING"]) / 2
        peaks = collect_peaks_higher(sfrm["IMG"])
        append!(df, DataFrame(:angles => Tuple(angles), :peak => peaks))
        println("$(basename(file)): collected $(length(peaks)) peaks")
    end
    select!(
        df,
        :angles =>
            ByRow(angles -> (angles[1]u"°", (angles[3]u"°", angles[2]u"°"))) => [:tth, :ϕω],
        :peak => [:coord, :size, :amplitude],
    )
end

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
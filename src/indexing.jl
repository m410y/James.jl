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

function estimate_reflex(
    wvec::AbstractVector,
    detector::Detector,
    gonio::Goniometer,
    detector_angle,
    gonio_angle,
    coord;
    cryst_pos::Point3 = Point3(0, 0, 0)u"m",
)
    p = rotate(detector, detector_angle...)(coord...)
    s = Vec3(wvec - norm(wvec) * normalize(p - cryst_pos))
    inv(gonio(gonio_angle...))(s)
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

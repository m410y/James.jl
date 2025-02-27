struct Frame{T,N}
    image::AbstractArray{T,N}
    setting::FrameSetting
end

struct Peak{N}
    intensity::Number
    baseline::Number
    coord::Vec{N}
    region::CartesianIndices{N}
    frame::Frame
end

function characterize_peak(frame::Frame, coord1, coord2)
    minmax_coord = extrema.(zip(coord1, coord2))
    min_coord = first.(minmax_coord)
    max_coord = last.(minmax_coord)
    region = CartesianIndex(min_coord...):CartesianIndex(max_coord...)
    baseline = median(frame.image[region])
    stat = sum(idx -> [1, idx.I...] * (frame.image[idx] - baseline), region)
    intensity = stat[begin]
    coord = stat[begin+1:end] / intensity
    Peak(intensity, baseline, Vec(coord...), region, frame)
end

function Base.show(io::IO, ::MIME"text/plain", peak::Peak{N}) where {N}
    print(io, "$N-dimentional Peak:\n")
    print(io, "  intensity: $(peak.intensity)\n")
    print(io, "  baseline: $(peak.baseline)\n")
    print(io, "  coordinate: $(Tuple(peak.coord))\n")
    print(io, "  region: $(peak.region.indices)")
end

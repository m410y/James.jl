struct Frame{T,N}
    image::AbstractArray{T,N}
    settings::FrameSettings
end

struct Peak{N}
    intensity::Number
    coord::Vec{N}
    indices::NTuple{N,AbstractRange}
    frame::Frame
end

function characterize_peak(frame::Frame, coord1, coord2)
    minmax_coord = extrema.(zip(coord1, coord2))
    min_coord = first.(minmax_coord)
    max_coord = last.(minmax_coord)
    region = CartesianIndex(min_coord...):CartesianIndex(max_coord...)
    stat = sum(idx -> [1, idx.I...] * image[idx] for idx in region)
    intensity = stat[begin]
    coord = stat[begin+1:end] / intensity
    Peak(intensity, coord, region.indices, frame)
end

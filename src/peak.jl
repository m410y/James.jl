struct Peak{N}
    int::Number
    coord::NTuple{N,Number}
    region::CartesianIndices{N}
end

function characterize_peak(image::AbstractArray, coord1, coord2)
    minmax_coord = extrema.(zip(coord1, coord2))
    min_coord = first.(minmax_coord)
    max_coord = last.(minmax_coord)
    region = CartesianIndex(min_coord...):CartesianIndex(max_coord...)
    stat = sum(idx -> [1, idx.I...] * image[idx] for idx in region)
    int = stat[begin]
    coord = stat[begin+1:end] / peak_int
    Peak(int, coord, region)
end

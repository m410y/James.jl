struct Detector{T}
    M::SMatrix{3,2,T,6}
    p::SVector{3,T}
end

const Coord{T<:Real} = SVector{2,T}

function coord(detector::Detector, line::BiVec)
    ex = direction(detector.M[:, 1])
    ey = direction(detector.M[:, 2])
    p0 = point(detector.p)
    x = (ey ∨ p0) ∨ line
    y = (ex ∨ p0) ∨ line
    t = (ex ∨ ey) ∨ line
    Coord(x / t, y / t)
end

point(detector::Detector, coord::Coord) = point(muladd(detector.M, coord, detector.p))

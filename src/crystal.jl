struct Crystal{T}
    UB::SMatrix{3,3,T,9}
    p::PVec
end

const Miller{T<:Real} = SVector{3,T}

miller(crystal::Crystal, q::PVec) = Miller(crystal.UB \ SA[q.e032, q.e013, q.e021])
direction(crystal::Crystal, hkl::Miller) = direction(crystal.UB * hkl)

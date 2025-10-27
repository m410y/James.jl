struct AngleCalculator{T<:Real}
    UB::Mat3{T}
    axes::NTuple{2,Vec3{T}}
    k::Vec3{T}
    function AngleCalculator{T}(
        UB::AbstractMatrix,
        axes,
        k::AbstractVector,
    ) where {T<:Real}
        new{T}(Mat3{T}(UB), Tuple(Vec3{T}(ax) for ax in axes), Vec3{T}(k))
    end
end
AngleCalculator(args...) = AngleCalculator{Float64}(args...)

function (acalc::AngleCalculator)(hkl::AbstractVector)
    q = acalc.UB * hkl
    try
        solve_equator_reflection(q, acalc.axes..., acalc.k)
    catch e
        if e isa DomainError 
            return nothing
        else
            rethrow(e)
        end
    end
end

function bragg_tth(acalc::AngleCalculator, hkl::AbstractVector)
    2asin(norm(acalc.UB * hkl) / 2norm(acalc.k))
end

default_hkl_iterator(acalc::AngleCalculator) = MillerIterator(acalc.UB'acalc.UB, 2norm(acalc.k))

function collect_angles(acalc::AngleCalculator{T}, hkls = default_hkl_iterator(acalc)) where {T<:Real}
    df = DataFrame("hkl" => Vec3i[], "2θ" => T[], "φ" => T[], "ω" => T[])
    for hkl in hkls
        angles = acalc(hkl)
        if isnothing(angles)
            continue
        end
        tth = bragg_tth(acalc, hkl)
        for (φ, ωs) in angles, ω in ωs
            push!(df, (hkl, tth, φ, ω))
        end
    end
    df
end

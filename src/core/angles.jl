angle(a::AbstractVector, b::AbstractVector) = atan(norm(cross(a, b)), dot(a, b))u"rad"

function axis_angle(b::AbstractVector, axis::AbstractVector, c::AbstractVector)
    b_a = b - axis * dot(axis, b)
    c_a = c - axis * dot(axis, c)
    atan(dot(axis, cross(b_a, c_a)), dot(b_a, c_a))u"rad"
end

function circle_intersection_angles(a, b, c)
    s = (a + b + c) / 2
    max_len, max_n = findmax((a, b, c))
    if s ≤ max_len
        if s < max_len
            return nothing
        end

        min_len, min_n = findmin((a, b, c))
        if min_len == 0
            angles = fill(pi / 2, 3)
            angles[min_n] = 0
        else
            angles = zeros(3)
            angles[max_n] = pi
        end

        α, β, γ = angles
    else
        sin_s = sin(s)
        sin_a = sin(s - a)
        sin_b = sin(s - b)
        sin_c = sin(s - c)

        α = 2atan(sqrt(sin_b * sin_c), sqrt(sin_s * sin_a))
        β = 2atan(sqrt(sin_c * sin_a), sqrt(sin_s * sin_b))
        γ = 2atan(sqrt(sin_a * sin_b), sqrt(sin_s * sin_c))
    end

    α * u"rad", β * u"rad", γ * u"rad"
end

function orient_angles(
    axis₁::AbstractVector,
    axis₂::AbstractVector,
    src::AbstractVector,
    dst::AbstractVector,
)
    SA₁ = angle(src, axis₁)
    A₁A₂ = angle(axis₁, axis₂)
    DA₂ = angle(dst, axis₂)

    ∠SA₁A₂ = axis_angle(src, axis₁, axis₂)
    ∠A₁A₂D = axis_angle(axis₁, axis₂, dst)
    Δα₁, Δα₂, _ = circle_intersection_angles(DA₂, SA₁, A₁A₂)

    [(∠SA₁A₂ + Δα₁, ∠A₁A₂D + Δα₂), (∠SA₁A₂ - Δα₁, ∠A₁A₂D - Δα₂)]
end

function reflection_angles(
    axis::AbstractVector,
    src::AbstractVector,
    dir::AbstractVector,
    θ::Number,
)
    SA = angle(src, axis)
    AD = angle(axis, dir)

    ∠SAD = axis_angle(src, axis, dir)
    Δα, _, _ = circle_intersection_angles(pi / 2 - θ, SA, AD)

    [∠SAD + Δα, ∠SAD - Δα]
end

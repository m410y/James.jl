struct Detector
    x::PVec
    y::PVec
    p::PVec
end

struct Coord{T<:Real}
    x::T
    y::T
end
Coord(x::Real, y::Real) = Coord(promote(x, y)...)
Coord(c::Coord) = c
Coord{T}(c::Coord) where {T<:Real} = Coord{T}(c.x, c.y)

promote_rule(::Type{Coord{T}}, ::Type{Coord{S}}) where {T<:Real,S<:Real} =
    Coord{promote_type(T,S)}

widen(::Type{Coord{T}}) where {T} = Coord{widen(T)}
float(::Type{Coord{T}}) where {T<:AbstractFloat} = Coord{T}
float(::Type{Coord{T}}) where {T} = Coord{float(T)}

isinteger(c::Coord) = isinteger(c.x) & isinteger(c.y)
isfinite(c::Coord) = isfinite(c.x) & isfinite(c.y)
isnan(c::Coord) = isnan(c.x) | isnan(c.y)
isinf(c::Coord) = isinf(c.x) | isinf(c.y)
iszero(c::Coord) = iszero(c.x) & iszero(c.y)

flipsign(c::Coord, y::Real) = ifelse(signbit(y), -c, c)
bswap(c::Coord) = Complex(bswap(c.x), bswap(c.y))

Base.:==(c1::Coord, c2::Coord) = (c1.x == c2.x) & (c1.y == c2.y)

Base.:+(c::Coord) = Coord(+c.x, +c.y)
Base.:-(c::Coord) = Coord(-c.x, -c.y)
Base.:+(c1::Coord, c2::Coord) = Coord(c1.x + c2.x, c1.y + c2.y)
Base.:-(c1::Coord, c2::Coord) = Coord(c1.x - c2.x, c1.y - c2.y)
Base.:*(s::Real, c::Coord) = Complex(s * c.x, s * c.y)
Base.:*(c::Coord, s::Real) = Complex(c.x * s, c.y * s)
Base.:/(c::Coord, s::Real) = Complex(c.x / s, c.y / s)

function line_coord(detector::Detector, line::BiVec)
    a = detector.p ∨ line
    b = detector.x ∨ detector.y ∨ line
    x = (a ∨ detector.y) / b
    y = (detector.x ∨ a) / b
    Coord(x, y)
end

function coord_point(detector::Detector, coord::Coord)
    detector.p + detector.x * coord.x + detector.y * coord.y
end

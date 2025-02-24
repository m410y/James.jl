function gauss(x::Number, x0::Number, σ::Number)
    exp(-((x - x0) / σ)^2 / 2) / sqrt(2π) / σ
end

function lorentz(x::Number, x0::Number, γ::Number)
    1 / (1 + ((x - x0) / γ)^2) / π / γ
end

# remove vectors translations, only points will be translated
(trans::Translation)(v::Vec3) = v
(aff::AffineMap)(v::Vec3) = aff.linear * v

module James

using LinearAlgebra
using Statistics
using Dates

using StaticArrays
using Unitful
using UnitfulEquivalences

using LeastSquaresOptim

using PGA3D
using SFRM
using P4P

const runit = u"keV"
const munit = u"mm"

include("detector.jl")
include("crystal.jl")
include("model.jl")
include("geometry.jl")

end # module James
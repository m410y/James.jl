module James

using LinearAlgebra

using GeometryBasics
using Rotations
using DataFrames

using Reexport
@reexport using JamesCore

export AngleCalculator, bragg_tth, collect_angles

include("angles.jl")

export D8Venture

include("d8venture.jl")

end

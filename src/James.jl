module James

using LinearAlgebra

using GeometryBasics
using DataFrames

using Reexport
@reexport using JamesCore

export AngleCalculator, bragg_tth, collect_angles

include("angles.jl")
include("propose.jl")

export D8Venture

include("d8venture.jl")

end

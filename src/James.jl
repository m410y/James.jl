module James

using LinearAlgebra, StaticArrays
using Unitful, UnitfulEquivalences
using GeometryBasics, CoordinateTransformations, Rotations

using Statistics
using IterTools

using DataFrames
using ImageFiltering, ImageMorphology
using FileIO

using Interpolations
using Integrals
using Optimization, FiniteDiff
using Zygote

include("core/utils.jl")
include("core/axis.jl")
include("core/angles.jl")
include("core/goniometer.jl")
include("core/detector.jl")
include("core/spectrum.jl")
include("core/crystal.jl")
include("core/model.jl")

include("peak.jl")
include("profile.jl")
include("refine.jl")

end

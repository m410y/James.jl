module James

using LinearAlgebra
using StaticArrays
using OffsetArrays
using Statistics
using Distributions
using IterTools
using Printf
using GeometryBasics
using CoordinateTransformations
using Rotations
using Unitful
using UnitfulEquivalences
using DataFrames
using ImageFiltering, ImageMorphology
using FileIO
using Interpolations
using Integrals
using Optimization, FiniteDiff, Zygote

include("core/utils.jl")
include("core/angles.jl")
include("core/axis.jl")
include("core/goniometer.jl")
include("core/detector.jl")
include("core/spectrum.jl")
include("core/sample.jl")

include("model.jl")
include("experiment.jl")
include("frame.jl")
include("load.jl")

include("refine/profile.jl")
include("refine/experiment.jl")

end

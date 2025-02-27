module James

using LinearAlgebra, StaticArrays
using OffsetArrays: Origin

using Statistics
using IterTools

using GeometryBasics, CoordinateTransformations, Rotations
using Unitful: Length, Energy
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
include("core/crystal.jl")
include("core/model.jl")

include("experiment.jl")
include("frame.jl")
include("load.jl")
include("profile.jl")
include("refine.jl")

end

module James

using LinearAlgebra, StaticArrays
using Unitful, UnitfulEquivalences
using GeometryBasics, CoordinateTransformations, Rotations

using Statistics
using DataFrames
using ImageFiltering, ImageMorphology
using FileIO

using Interpolations
using Integrals
using Optimization

# using Crystalline
# using Statistics
# using DataFrames
# using Images
# using LsqFit
# using FFTW

include("utils.jl")
include("axis.jl")
include("angles.jl")
include("goniometer.jl")
include("detector.jl")
include("spectrum.jl")
include("profile.jl")
include("crystal.jl")
include("model.jl")
include("refine.jl")
include("indexing.jl")

# include("process.jl")
# include("predict_reflections.jl")
# include("peak_profile.jl")

end # module James

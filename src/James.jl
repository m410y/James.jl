module James

using LinearAlgebra

using GeometryBasics
using Rotations
using DataFrames
using Unitful

using Statistics
import ForwardDiff
using Optimization: OptimizationFunction, OptimizationProblem, solve, AutoForwardDiff
using OptimizationOptimJL: Newton
using DifferentiationInterface: SecondOrder

using Reexport
@reexport using JamesCore

export AngleCalculator, bragg_tth, collect_angles

include("angles.jl")

export D8Venture

include("d8venture.jl")

export fit_peak

include("fit.jl")

end

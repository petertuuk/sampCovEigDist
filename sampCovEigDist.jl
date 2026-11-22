# Kept so that existing scripts using `include("sampCovEigDist.jl")` still work.
# The code now lives in src/SampCovEigDist.jl as a module.
include(joinpath(@__DIR__, "src", "SampCovEigDist.jl"))
using .SampCovEigDist

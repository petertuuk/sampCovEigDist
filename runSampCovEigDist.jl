## Setup
v = VERSION
if v >= v"0.7"
    using LinearAlgebra
end
include("sampCovEigDist.jl")

import DelimitedFiles, Plots, BenchmarkTools

# Calculate Test Case
c = .1
evweight = [.2 1; .4 3; .4 10];
density, zvec, cdf  = sampCovEigDist(evweight, c)

# Compare Result to Exepected Result
expected_density = DelimitedFiles.readdlm("expected_density.csv", ',', Any, '\n')
Plots.plot(expected_density[:,1], expected_density[:,2],linewidth=4,label="Expected Density")
Plots.plot!(zvec, density,ylims=(0,.5), xlims=(0,15),linewidth=2,label="Calculated Density")
Plots.xlabel!("Eigenvalue Magnitude")
Plots.ylabel!("Eigenvalue Distribution Density")
Plots.gui()

BenchmarkTools.@btime sampCovEigDist(evweight, c)

## profile
# v = VERSION
# if v >= v"0.7"
#     using Profile
# end
# Profile.init()
# Profile.clear()
# @profile sampCovEigDist(evweight, c)
# f = open("profile.prof", "w")
# Profile.print(f)
# close(f)
# Profile.clear_malloc_data()
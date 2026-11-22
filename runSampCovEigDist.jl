## Setup
include("sampCovEigDist.jl")

import DelimitedFiles, Plots, BenchmarkTools

# Calculate test case
c = 0.1
evweight = [.2 1; .4 3; .4 10]   # columns: [weight eigenvalue]
density, zvec, cdf = sampCovEigDist(evweight, c)

# Compare result to expected result
expected_density = DelimitedFiles.readdlm("expected_density.csv", ',', Float64)
Plots.plot(expected_density[:, 1], expected_density[:, 2], linewidth=4, label="Expected Density")
Plots.plot!(zvec, density, ylims=(0, .5), xlims=(0, 15), linewidth=2, label="Calculated Density")
Plots.xlabel!("Eigenvalue Magnitude")
Plots.ylabel!("Eigenvalue Distribution Density")
Plots.gui()

BenchmarkTools.@btime sampCovEigDist($evweight, $c)

## Profile (uncomment to use)
# using Profile
# Profile.clear()
# @profile sampCovEigDist(evweight, c)
# open("profile.prof", "w") do f
#     Profile.print(f)
# end

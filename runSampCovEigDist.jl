## Setup
v = VERSION
if v >= v"0.7"
    using LinearAlgebra
end
include("sampCovEigDist.jl")
if v >= v"0.7"
    ee = range(5, stop=75, length=400)
else
    ee = linspace(5, 75, 400)
end

evweight = [2*length(ee) 1; [fill(1.0,length(ee)) ee]]
const c = .5

## Calc/Plot
doPlot = true
if doPlot
    density1, zvec, cdf = sampCovEigDist(evweight, c)
    evweight2 = [.8 1; .1 4; .1 10];
    density2, zvec2, cdf2  = sampCovEigDist(evweight2, c, zvec)
    evweight3 = [1 1]
    density3, zvec3, cdf3  = sampCovEigDist(evweight3, c, zvec)

    using Plots;pyplot();
    plot( zvec, density1 .+ 1e-16)
    plot!(zvec, density2 .+ 1e-16)
    plot!(zvec, density3 .+ 1e-16, ylims=(1e-3,10), xlims=(.07,100), xscale=:log10, yscale=:log10)
    gui()
end

sampCovEigDist(evweight, c)

using BenchmarkTools
@btime sampCovEigDist(evweight, c)

## profile
v = VERSION
if v >= v"0.7"
    using Profile
end
Profile.init()
Profile.clear()
@profile sampCovEigDist(evweight, c)
f = open("profile.prof", "w")
Profile.print(f)
close(f)
# using ProfileView
# ProfileView.view()

Profile.clear_malloc_data()
sampCovEigDist(evweight, c)

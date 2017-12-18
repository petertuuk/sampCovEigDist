## Setup
include("sampCovEigDist.jl")

ee = linspace(5,75,400)
evweight = [2*length(ee) 1;[fill(1.0,length(ee)) ee]]
const c = .5

## Calc/Plot
doPlot = true
if doPlot
    density, zvec, cdf = sampCovEigDist(evweight,c)
    evweight2 = [.8 1;.1 4;.1 10];
    density2, zvec2, cdf2  = sampCovEigDist(evweight2,c,zvec)
    density3, zvec3, cdf3  = sampCovEigDist([1 1],c,zvec)
    using PyPlot
    doLog = true
    if doLog
      semilogy(zvec,density)
      semilogy(zvec,density2)
      semilogy(zvec,density3)
      ylim(1e-3, 10)
    else
      plot(zvec,density)
      plot(zvec,density2)
      plot(zvec,density3)
    end
    show()
    # Time
end

sampCovEigDist(evweight,c)

@time sampCovEigDist(evweight,c)

## profile
Profile.init()
Profile.clear_malloc_data
Profile.clear()
@profile sampCovEigDist(evweight,c)
f = open("profile.prof","w")
Profile.print(f)
close(f)
# using ProfileView
# ProfileView.view()

function sampCovEigDist(evweight, c, zvec2in=0)
# (density, zvec2, cdf) = getDensity(evweight, c, zvec2in=0)
# This function computes the expected sample spectral distribution
#
# Inputs:
#   - evweight: a two-column matrix in which the first column is (samples of) the true eigenvalue distribution
#     and the second column is the mass at each corresponding eigenvalue.
#   - c: the sampling ratio (number of samples to dimensionality of the random variable)
#   - zvec2in (optional): the points at which the return distribution will be discretized
# Outputs:
#   - density: the eigenvalue density
#   - zvec2: the eigenvalue domain at which the density is specified
#   - cdf: a cumulative sum of the density
#
# Uses Stieltjes transform convergence and inversion
# Ref: Bai and Silverstein, "No Eigenvalues Outside the Support of the Limiting
#      Spectral Distribution of Large Dimensional Sample Covariance Matrices"
#      Ann. Prob., 26 (1), 316-345.
#      https://projecteuclid.org/euclid.aop/1022855421

  v = VERSION

  nMbar = 5000 # sufficient for most
  nZvec = 5000

  evweight[:,1] = evweight[:,1]./sum(evweight[:,1])
  if v.major>0 || v.minor>6
    evweight = round.(evweight,digits=6)
  else
    evweight = round.(evweight,6)
  end

  #generate functions to compute common integrals
  cumputeInteg,cumputeIntegSq = makeIntegralFuncs(evweight)

  # find start and end-points where intervals where eigenvalues live
  intervalStarts,intervalEnds = getIntervals(nMbar,evweight,cumputeInteg)

  zvec2::Vector{Float64} = []
  if length(zvec2in)>1
    zvec2 = zvec2in;
    nZvec = length(zvec2)
  else
    ll = 3;

    if v >= v"0.7"
      zvec2 = range(0,stop=(1.01*intervalEnds[end]).^(1/ll),length=nZvec+1).^ll
    else
      zvec2 = linspace(0,(1.01*intervalEnds[end]).^(1/ll),nZvec+1).^ll
    end
    zvec2 = zvec2[2:end];
  end

  density = fill(0.0,nZvec);
  wentBadly = false
  guess0 = -5.0 + 10.0im
  guess = 0.0im
  fminres = 0.0im
  fminresOld = 0.0im
  cfact = 1/min(c,1)/pi

  # Calculate which points in eigenvalue sample space are in the support
  islt = fill(0.0,size(zvec2))
  ielt = fill(0.0,size(zvec2))
  intervalStarts99::Vector{Float64} = intervalStarts*.99
  intervalEnds101::Vector{Float64} = intervalEnds*1.01
  for iz in eachindex(zvec2)
    for ii1 in eachindex(intervalStarts99)
      if zvec2[iz]>intervalStarts99[ii1]
        islt[iz] += 1
      end
      if zvec2[iz]>intervalEnds101[ii1]
        ielt[iz] += 1
      end
    end
  end

  #Main Loop: loop over the domain
  zast::Float64 = 0
  for iz in eachindex(zvec2)
    zast = zvec2[iz]
    if islt[iz] > ielt[iz]

      # Define the objective function and its derivative
      @inline function mobj(m::Complex{Float64})
        m + 1/( zast - cumputeInteg(m) )
      end
      @inline function mobjDiff(m::Complex{Float64})
        1 - cumputeIntegSq(m) / ( zast - cumputeInteg(m) )^2
      end

      # Make the guess using linear prediction from last two results
      if !wentBadly && imag(fminres) > 1e-8
        guess = fminres + fminres - fminresOld
      else
        guess = guess0
      end

      # Perform the optimization to find the density at this point
      fminresOld = fminres
      fminres = newton(mobj, mobjDiff, guess);
      density[iz] = abs(imag(fminres)) * cfact

      # Check the result to ensure consistency
      z_recalc = -1/fminres + cumputeInteg(fminres);
      z_resid = abs(z_recalc - zast);
      if z_resid > 1e-3 && z_resid/abs(zast) > 1e-3
        wentBadly = true;
        print("some divergence in z: $(z_resid/abs(zast))");
      else
        wentBadly = false
      end
    end
  end

  # Calculate CDF using trapezoidal integration
  cdf = trapz(zvec2,density)
  if cdf[end]>1.02 || cdf[end]<1/1.02
    warn("cdf was $(cdf[end])!\n")
  end
  cdf = cdf/cdf[end]

  return density, zvec2, cdf
end

@inline function newton(f::Function, df::Function, x0::Complex{Float64}, tol=1e-6, kmax=100)
# Textbook Newton solver
  x::Complex{Float64} = 0.0
  ex = tol + 1.0
  k::Int16 = 0;
  while ex > tol && k <= kmax
    k += 1
    x  = x0 - f(x0) / df(x0)
    ex = abs(x - x0) / abs(x)
    x0 = x
  end
  return x
end

function trapz(x,y)
# Trapezoid integration
  area = fill(0.0,size(x))
  lx = length(x)
  for i in 1:lx-1
    area[i+1] = area[i] + (x[i+1] - x[i]) * (y[i+1] + y[i]) * .5
  end
  return area
end

function makeMbarVec(evweight,nMbar)
  startm = min(log10(1/maximum(evweight[:,2])/10),-3)
  endm = max(-log10(1/minimum(evweight[:,2])),3)
  if v.major>0 || v.minor>6
    mbarvecHalf = 10 .^ range(startm,stop=endm,length=Int(nMbar/2))
  else
    mbarvecHalf = logspace(startm,endm,Int(nMbar/2))
  end
  mbarvec = sort([-mbarvecHalf;mbarvecHalf])
end

function makeIntegralFuncs(evweight)
# Generate the functions that will compute commonly used integral and
# its partial derivative
  ut::Array{Float64,1} = unique(evweight[:,2]);
  ht::Array{Float64,1} = fill(0.0,size(ut))
  for iu in eachindex(ut)
    ixEv = evweight[:,2].== ut[iu]
    ht[iu] = sum(evweight[ixEv,1])
  end

  uh::Array{Float64,1} = ut.*ht;
  uhSq::Array{Float64,1} = (ut.^2) .* ht;
  # uhSq::Array{Float64,1} = ut .* ut .* ht;
  @inline function cumputeInteg(mbar)
    funcsum::Complex{Float64} = 0
    @inbounds for ii2 in eachindex(ut)
      funcsum += uh[ii2]/(1+ut[ii2]*mbar);
    end
    funcsum *= c
    return funcsum
  end
  @inline function cumputeIntegSq(mbar)
    funcsum::Complex{Float64} = 0
    @inbounds for ii3 in eachindex(ut)
      funcsum += uhSq[ii3]/(1+ut[ii3]*mbar)^2;
    end
    funcsum *= c
    return funcsum
  end

  return cumputeInteg,cumputeIntegSq
end

function getIntervals(nMbar,evweight,cumputeInteg)
# find start and stop values for each interval in which
# the spectral distribution will have non-zero mass
  mbarvec = makeMbarVec(evweight,nMbar)

  zvec = fill(0.0,size(mbarvec))
  @simd for i1 in eachindex(mbarvec)
    mbar = mbarvec[i1];
    zvec[i1] = -1/mbar + cumputeInteg(mbar);
  end

  eigInds::Array{Int32,1} = fill(0.0,size(evweight,1))
  evDiff = 0.0
  for uu in 1:size(evweight,1)
    evInv::Float64 = -1/evweight[uu,2]
    minval = 9999.9
    minind = 0
    for iu = eachindex(mbarvec)
      evDiff = abs(evInv-mbarvec[iu])
      if  evDiff < minval
        minval = abs(evInv-mbarvec[iu])
        minind = iu
      end
    end
    eigInds[uu] = minind
  end
  isPeak = unique(eigInds);

  diffzvec = diff(zvec)
  picks = diffzvec .> 0
  picks[isPeak.-1] .= false
  picks[isPeak] .= false
  picks[isPeak.+1] .= false

  if v >= v"0.7"
    findPicks = findall(picks)
    dfpUp = findall(diff([-Inf;findPicks]).>1)
    dfpDown = findall(diff([findPicks;Inf]).>1)
  else
    findPicks = find(picks)
    dfpUp = find(diff([-Inf;findPicks]).>1)
    dfpDown = find(diff([findPicks;Inf]).>1)
  end

  zvs = zvec[findPicks[dfpUp]]
  zve = zvec[findPicks[dfpDown]]
  if v >= v"0.7"
    complementOfSupport = sort([zvs zve],dims=1)
  else
    complementOfSupport = sort([zvs zve],1)
  end
  intervalStarts = complementOfSupport[1:end-1,2];
  intervalEnds = complementOfSupport[2:end,1];
  if length(intervalStarts)>0
    intLength = intervalEnds-intervalStarts
    intBreak = intervalStarts[2:end]-intervalEnds[1:end-1]
    if length(intBreak)>0
      intervalStarts = intervalStarts[[true;intBreak.>0]]
      intervalEnds = intervalEnds[[intBreak.>0;true]]
    end
  end

  intervalStarts[intervalStarts.<0] .= 0;

  if any( intervalStarts[2:end] .<= intervalEnds[1:end-1] )
    dse = intervalStarts[2:end] - intervalEnds[1:end-1]
    takeStarts = [1;find(dse.>0)+1];
    takeEnds = [find(dse.>0);length(intervalStarts)];

    intervalStarts = intervalStarts[takeStarts];
    intervalEnds = intervalEnds[takeEnds];
  end
  return intervalStarts,intervalEnds
end

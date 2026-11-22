module SampCovEigDist

export sampCovEigDist

"""
    density, zvec, cdf = sampCovEigDist(evweight, c, zvec_in=nothing; nMbar=5000, nZvec=5000)

Compute the limiting (large-dimension) eigenvalue density of a sample covariance
matrix, given the eigenvalue distribution of the true covariance matrix.

# Inputs
- `evweight`: two-column matrix describing the true (population) eigenvalue
  distribution as a set of point masses. Column 1 is the weight (mass) of each
  eigenvalue, column 2 is the eigenvalue itself. Weights are normalized to sum to
  one; repeated eigenvalues are merged. The input is not modified.
- `c`: dimension-to-sample ratio, `c = p/n` (dimension of the random vector
  divided by the number of samples).
- `zvec_in` (optional): grid of eigenvalues at which to evaluate the density.
  If omitted (or `0`/`nothing`), a grid of `nZvec` points covering the support
  is generated.

# Keywords
- `nMbar`: number of points used to trace out the support.
- `nZvec`: number of points in the generated output grid.

# Outputs
- `density`: the eigenvalue density at each point of `zvec`
- `zvec`: the grid of eigenvalues
- `cdf`: cumulative distribution, by trapezoidal integration of `density`,
  normalized to end at one

For `c > 1` the sample covariance matrix has a point mass of `1 - 1/c` at zero.
That mass is not included: the returned density and cdf describe the
distribution of the nonzero eigenvalues.

Uses Stieltjes transform convergence and inversion.
Ref: Bai and Silverstein, "No Eigenvalues Outside the Support of the Limiting
Spectral Distribution of Large Dimensional Sample Covariance Matrices",
Ann. Prob., 26 (1), 316-345. https://projecteuclid.org/euclid.aop/1022855421
"""
function sampCovEigDist(evweight::AbstractMatrix{<:Real}, c::Real, zvec_in=nothing;
                        nMbar::Integer=5000, nZvec::Integer=5000)
    size(evweight, 2) == 2 || throw(ArgumentError("evweight must have two columns: [weight eigenvalue]"))
    c > 0 || throw(ArgumentError("c must be positive"))
    all(evweight[:, 2] .> 0) || throw(ArgumentError("eigenvalues (column 2) must be positive"))
    all(evweight[:, 1] .>= 0) && sum(evweight[:, 1]) > 0 ||
        throw(ArgumentError("weights (column 1) must be nonnegative and not all zero"))

    spec = PopulationSpectrum(evweight, c)

    # find start and end points of the intervals where the eigenvalues live
    intervalStarts, intervalEnds = getIntervals(spec, nMbar)

    # grid on which to evaluate the density
    if zvec_in === nothing || (zvec_in isa Number)
        ll = 3   # cube-root spacing puts more points near zero
        zvec = collect(range(0, stop=(1.01 * intervalEnds[end])^(1 / ll), length=nZvec + 1) .^ ll)[2:end]
    else
        zvec = collect(Float64, zvec_in)
    end

    density = zeros(length(zvec))
    cfact = 1 / min(c, 1) / pi

    # which grid points are (within a small margin of) the support
    starts99 = intervalStarts .* 0.99
    ends101 = intervalEnds .* 1.01
    inSupport = [count(<(z), starts99) > count(<(z), ends101) for z in zvec]

    # Main loop: solve the Silverstein equation for the Stieltjes transform at each z
    guess0 = complex(-5.0, 10.0)
    fminres = zero(ComplexF64)
    fminresOld = zero(ComplexF64)
    wentBadly = false
    nBad = 0
    worstResid = 0.0
    for iz in eachindex(zvec)
        inSupport[iz] || continue
        z = zvec[iz]

        # guess using linear prediction from the last two results
        guess = (!wentBadly && imag(fminres) > 1e-8) ? 2fminres - fminresOld : guess0

        fminresOld = fminres
        fminres = solveStieltjes(spec, z, guess)
        density[iz] = abs(imag(fminres)) * cfact

        # check the result by mapping back to z
        zResid = abs(-1 / fminres + integ(spec, fminres) - z)
        if zResid > 1e-3 && zResid / abs(z) > 1e-3
            wentBadly = true
            nBad += 1
            worstResid = max(worstResid, zResid / abs(z))
        else
            wentBadly = false
        end
    end
    if nBad > 0
        @warn "Stieltjes solver did not converge at $nBad of $(count(inSupport)) points (worst relative z residual $worstResid)"
    end

    cdf = trapz(zvec, density)
    if !(1 / 1.02 <= cdf[end] <= 1.02)
        @warn "Density integrates to $(cdf[end]) rather than 1; the cdf has been renormalized. Try a finer zvec."
    end
    cdf ./= cdf[end]

    return density, zvec, cdf
end

"""
Population eigenvalue distribution, merged into unique eigenvalues `t` with
masses `h`, along with the precomputed products used in the integrals.
"""
struct PopulationSpectrum
    t::Vector{Float64}
    th::Vector{Float64}     # t .* h
    t2h::Vector{Float64}    # t.^2 .* h
    c::Float64
end

function PopulationSpectrum(evweight::AbstractMatrix{<:Real}, c::Real)
    w = Float64.(evweight[:, 1])
    ev = Float64.(evweight[:, 2])
    w ./= sum(w)
    t = sort(unique(ev[w .> 0]))
    h = [sum(w[ev .== ti]) for ti in t]
    return PopulationSpectrum(t, t .* h, t .^ 2 .* h, Float64(c))
end

"c * ∫ t / (1 + t m) dH(t)"
@inline function integ(s::PopulationSpectrum, m)
    acc = zero(complex(float(typeof(m))))
    @inbounds for i in eachindex(s.t)
        acc += s.th[i] / (1 + s.t[i] * m)
    end
    return s.c * acc
end

"Both c ∫ t/(1+tm) dH and c ∫ t²/(1+tm)² dH in one pass."
@inline function integPair(s::PopulationSpectrum, m::ComplexF64)
    a = zero(ComplexF64)
    b = zero(ComplexF64)
    @inbounds for i in eachindex(s.t)
        r = 1 / (1 + s.t[i] * m)
        a += s.th[i] * r
        b += s.t2h[i] * r * r
    end
    return s.c * a, s.c * b
end

"""
Newton's method on f(m) = m + 1/(z - c∫t/(1+tm)dH), the Silverstein equation
for the companion Stieltjes transform m at the real point z.
"""
function solveStieltjes(s::PopulationSpectrum, z::Float64, m0::ComplexF64; tol=1e-6, kmax=100)
    m = m0
    for _ in 1:kmax
        I1, I2 = integPair(s, m)
        d = z - I1
        f = m + 1 / d
        df = 1 - I2 / d^2
        mNew = m - f / df
        converged = abs(mNew - m) <= tol * abs(mNew)
        m = mNew
        converged && break
    end
    return m
end

"Cumulative trapezoidal integration of y over x."
function trapz(x, y)
    area = zeros(length(x))
    for i in 1:length(x)-1
        area[i+1] = area[i] + (x[i+1] - x[i]) * (y[i+1] + y[i]) / 2
    end
    return area
end

"Log-spaced grid of real (negative and positive) values of m used to trace z(m)."
function makeMbarVec(t::Vector{Float64}, nMbar::Integer)
    startm = min(log10(1 / maximum(t) / 10), -3)
    endm = max(log10(minimum(t)), 3)
    half = 10 .^ range(startm, stop=endm, length=nMbar ÷ 2)
    return sort([-half; half])
end

"Index of the element of sorted vector `v` nearest to `x`."
function nearestIndex(v::AbstractVector, x)
    i = searchsortedfirst(v, x)
    i > length(v) && return length(v)
    i == 1 && return 1
    return (x - v[i-1]) <= (v[i] - x) ? i - 1 : i
end

"""
Find the start and end of each interval on which the limiting spectral
distribution has nonzero mass. Uses the Silverstein & Choi characterization:
outside the support, z(m) = -1/m + c∫t/(1+tm)dH is increasing in real m.
"""
function getIntervals(spec::PopulationSpectrum, nMbar::Integer)
    mbarvec = makeMbarVec(spec.t, nMbar)
    zvec = [-1 / m + real(integ(spec, m)) for m in mbarvec]

    # z(m) has poles at m = -1/t; exclude grid points next to them
    isPeak = unique([nearestIndex(mbarvec, -1 / ti) for ti in spec.t])

    picks = diff(zvec) .> 0
    np = length(picks)
    for k in isPeak, j in (k - 1, k, k + 1)
        1 <= j <= np && (picks[j] = false)
    end

    findPicks = findall(picks)
    isempty(findPicks) && error("Could not locate the support of the spectral distribution; try increasing nMbar")
    dfpUp = findall(diff([-Inf; findPicks]) .> 1)
    dfpDown = findall(diff([findPicks; Inf]) .> 1)

    # each run of increasing z(m) is an interval in the complement of the support
    complementOfSupport = sort([zvec[findPicks[dfpUp]] zvec[findPicks[dfpDown]]], dims=1)
    intervalStarts = complementOfSupport[1:end-1, 2]
    intervalEnds = complementOfSupport[2:end, 1]
    isempty(intervalStarts) && error("Could not locate the support of the spectral distribution; try increasing nMbar")

    intBreak = intervalStarts[2:end] - intervalEnds[1:end-1]
    if !isempty(intBreak)
        intervalStarts = intervalStarts[[true; intBreak .> 0]]
        intervalEnds = intervalEnds[[intBreak .> 0; true]]
    end

    intervalStarts[intervalStarts .< 0] .= 0

    # merge any intervals that overlap
    if any(intervalStarts[2:end] .<= intervalEnds[1:end-1])
        dse = intervalStarts[2:end] - intervalEnds[1:end-1]
        takeStarts = [1; findall(dse .> 0) .+ 1]
        takeEnds = [findall(dse .> 0); length(intervalStarts)]
        intervalStarts = intervalStarts[takeStarts]
        intervalEnds = intervalEnds[takeEnds]
    end
    return intervalStarts, intervalEnds
end

end # module

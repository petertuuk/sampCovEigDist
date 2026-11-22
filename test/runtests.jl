using Test, LinearAlgebra, Random, Statistics, DelimitedFiles

if !isdefined(Main, :SampCovEigDist)
    include(joinpath(@__DIR__, "..", "src", "SampCovEigDist.jl"))
end
using .SampCovEigDist

# Marchenko-Pastur density (true covariance = identity), c = p/n <= 1
function mpDensity(x, c)
    a, b = (1 - sqrt(c))^2, (1 + sqrt(c))^2
    return a < x < b ? sqrt((b - x) * (x - a)) / (2pi * c * x) : 0.0
end

# linear interpolation of (xs, ys) at x, zero outside
function interp(x, xs, ys)
    i = searchsortedlast(xs, x)
    (i < 1 || i >= length(xs)) && return 0.0
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    return ys[i] * (1 - t) + ys[i+1] * t
end

# empirical eigenvalues of sample covariance matrices with true eigenvalues ev
function simulatedEigs(ev, n; trials=20, seed=1)
    rng = MersenneTwister(seed)
    p = length(ev)
    out = Float64[]
    for _ in 1:trials
        X = sqrt.(ev) .* randn(rng, p, n)
        append!(out, eigvals(Symmetric(X * X' / n)))
    end
    return out
end

@testset "SampCovEigDist" begin

    @testset "regression: expected_density.csv" begin
        expected = readdlm(joinpath(@__DIR__, "..", "expected_density.csv"), ',', Float64)
        density, zvec, cdf = sampCovEigDist([.2 1; .4 3; .4 10], 0.1)
        calc = [interp(x, zvec, density) for x in expected[:, 1]]
        @test maximum(abs.(calc .- expected[:, 2])) < 1e-6
        @test cdf[end] ≈ 1
    end

    @testset "Marchenko-Pastur, c = $c" for c in (0.05, 0.1, 0.3, 0.5, 0.9)
        density, zvec, cdf = sampCovEigDist([1.0 1.0], c)
        # compare away from the square-root edges, where the grid limits accuracy
        a, b = (1 - sqrt(c))^2, (1 + sqrt(c))^2
        inner = (zvec .> a + 0.02 * (b - a)) .& (zvec .< b - 0.02 * (b - a))
        @test maximum(abs.(density[inner] .- mpDensity.(zvec[inner], c))) < 1e-4
        # (almost) no mass reported well outside the support
        @test all(density[(zvec .< 0.9a) .| (zvec .> 1.1b)] .< 1e-6)
        # density integrates to one without needing renormalization
        @test isapprox(sum(diff(zvec) .* (density[1:end-1] .+ density[2:end]) ./ 2), 1, atol=0.01)
    end

    @testset "c > 1 gives the distribution of the nonzero eigenvalues" begin
        c = 2.0
        density, zvec, cdf = sampCovEigDist([1.0 1.0], c)
        # nonzero eigenvalues of a p x p (p > n) sample covariance are the eigenvalues of
        # the n x n companion matrix, which is Marchenko-Pastur with ratio 1/c, scaled by c
        ref(x) = mpDensity(x / c, 1 / c) / c
        a, b = (1 - sqrt(c))^2, (1 + sqrt(c))^2
        inner = (zvec .> a + 0.02 * (b - a)) .& (zvec .< b - 0.02 * (b - a))
        @test maximum(abs.(density[inner] .- ref.(zvec[inner]))) < 1e-4
    end

    @testset "agrees with Monte Carlo" begin
        cases = [([.8 1; .1 4; .1 10], 400, 800),       # c = 0.5
                 ([.2 1; .4 3; .4 10], 200, 2000)]      # c = 0.1
        for (ew, p, n) in cases
            ev = vcat([fill(ew[k, 2], round(Int, ew[k, 1] * p)) for k in 1:size(ew, 1)]...)
            sims = simulatedEigs(ev, n)
            density, zvec, cdf = sampCovEigDist(ew, p / n)
            gap = maximum(abs(mean(sims .<= zvec[i]) - cdf[i]) for i in 1:25:length(zvec))
            @test gap < 0.01
        end
    end

    @testset "uses the c that is passed in" begin
        # regression for an old bug where a global `c` was used instead of the argument
        ew = [.2 1; .4 3; .4 10]
        d1, z1, _ = sampCovEigDist(ew, 0.1)
        d5, z5, _ = sampCovEigDist(ew, 0.5)
        @test z1[end] < z5[end]   # wider support for larger c
        @test (@test_logs min_level = Base.CoreLogging.Warn sampCovEigDist(ew, 0.5)) isa Tuple
    end

    @testset "input handling" begin
        ew = [2.0 1; 2.0 4]
        ewCopy = copy(ew)
        sampCovEigDist(ew, 0.3)
        @test ew == ewCopy                           # input not modified

        dInt, zInt, _ = sampCovEigDist([1 1; 1 4], 0.3)   # integer input works
        dF, zF, _ = sampCovEigDist([0.5 1; 0.5 4], 0.3)
        @test dInt ≈ dF && zInt ≈ zF

        # repeated eigenvalues are merged
        dRep, _, _ = sampCovEigDist([0.25 1; 0.25 1; 0.5 4], 0.3)
        @test dRep ≈ dF

        # user-supplied grid
        grid = range(0.01, 10, length=300)
        dG, zG, _ = sampCovEigDist([0.5 1; 0.5 4], 0.3, grid)
        @test zG == collect(grid) && length(dG) == 300

        @test_throws ArgumentError sampCovEigDist([1.0 1.0 1.0], 0.5)
        @test_throws ArgumentError sampCovEigDist([1.0 -1.0], 0.5)
        @test_throws ArgumentError sampCovEigDist([1.0 1.0], 0)
    end

    @testset "many-eigenvalue example (overlapping intervals)" begin
        ee = range(5, stop=75, length=400)
        ew = [2 * length(ee) 1; [fill(1.0, length(ee)) ee]]
        density, zvec, cdf = sampCovEigDist(ew, 0.5)
        @test all(isfinite, density) && cdf[end] ≈ 1
    end
end

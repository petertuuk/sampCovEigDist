Calculate the distribution of the eigenvalues of sample covariance matrices

Given the eigenvalues of a true covariance matrix and the ratio `c = p/n`
(dimension over number of samples), computes the limiting density of the
eigenvalues of the sample covariance matrix.

Uses Stieltjes transform convergence and inversion
Ref: Bai and Silverstein, "No Eigenvalues Outside the Support of the Limiting
Spectral Distribution of Large Dimensional Sample Covariance Matrices"
Ann. Prob., 26 (1), 316-345.
https://projecteuclid.org/euclid.aop/1022855421

[![DOI](https://zenodo.org/badge/91010751.svg)](https://zenodo.org/badge/latestdoi/91010751)

## Usage

Requires Julia 1.6 or later.

```julia
include("sampCovEigDist.jl")          # or: using SampCovEigDist (after Pkg.develop)

evweight = [.2 1; .4 3; .4 10]       # columns: [weight eigenvalue]
c = 0.1                              # p/n
density, zvec, cdf = sampCovEigDist(evweight, c)

# evaluate on your own grid
density, zvec, cdf = sampCovEigDist(evweight, c, range(0.01, 15, length=1000))
```

For `c > 1` the point mass of `1 - 1/c` at zero is left out, and the result
describes the nonzero eigenvalues.

`runSampCovEigDist.jl` plots the example above against `expected_density.csv`
(needs Plots and BenchmarkTools).

## Tests

```
julia test/runtests.jl
```

or, from the package REPL in this folder, `activate .` then `test`. The tests
compare against `expected_density.csv`, the closed-form Marchenko-Pastur
density, and Monte Carlo simulations of sample covariance matrices.

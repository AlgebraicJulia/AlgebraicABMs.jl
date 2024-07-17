# AlgebraicABMs.jl benchmarks

This directory contains benchmarks for different parts of AlgebraicABMs. To run all the
benchmarks, launch `julia --project=benchmark` and enter:

``` julia
using PkgBenchmark
import AlgebraicABMs

benchmarkpkg(AlgebraicABMs)
```

To run a specific set of benchmarks, use the `script` keyword argument, for
example:

``` julia
benchmarkpkg(AlgebraicABMs; script="benchmark/ABMs.jl")
```

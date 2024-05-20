# AlgebraicABMs.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://AlgebraicJulia.github.io/AlgebraicABMs.jl/stable)
[![Development Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://AlgebraicJulia.github.io/AlgebraicABMs.jl/dev)
[![Code Coverage](https://codecov.io/gh/AlgebraicJulia/AlgebraicABMs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AlgebraicJulia/AlgebraicABMs.jl)
[![CI/CD](https://github.com/AlgebraicJulia/AlgebraicABMs.jl/actions/workflows/julia_ci.yml/badge.svg)](https://github.com/AlgebraicJulia/AlgebraicABMs.jl/actions/workflows/julia_ci.yml)

Important example files: `docs/literate/game_of_life.jl` and `test/ABMs.jl`

## Caveats

[Fleck.jl](https://github.com/adolgert/Fleck.jl) is not officially released, so you must clone that repo and `dev` it.

This will need to be addressed before the tests can be run / documentation can be built automatically by Github actions.

## 🛠️ Usage

To locally build the documentation and the literate code examples, run the following in the command line:
```
julia --project=docs -e "using AlgebraicABMs, LiveServer; servedocs(literate_dir=\"docs/literate\",skip_dir=\"docs/src/generated\")"
```

To locally run the test suite, run the following command
```
julia --project=test test/runtests.jl
```

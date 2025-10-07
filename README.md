# AlgebraicABMs.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://AlgebraicJulia.github.io/AlgebraicABMs.jl/stable)
[![Development Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://AlgebraicJulia.github.io/AlgebraicABMs.jl/dev)
[![Code Coverage](https://codecov.io/gh/AlgebraicJulia/AlgebraicABMs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AlgebraicJulia/AlgebraicABMs.jl)
[![CI/CD](https://github.com/AlgebraicJulia/AlgebraicABMs.jl/actions/workflows/julia_ci.yml/badge.svg)](https://github.com/AlgebraicJulia/AlgebraicABMs.jl/actions/workflows/julia_ci.yml)

Important example files: `docs/literate/game_of_life.jl` and `test/ABMs.jl`

## 🛠️ Usage

To install the dependencies and register this package so that it can be imported, run the following command line:

```
julia -e 'using Pkg; Pkg.develop(path=".")'
```

but it might be better to use `--project` at this time rather than registering the project globally.

To locally build the documentation and the literate code examples, first register this package and install the dependencies:

```
julia --project=docs -e 'using Pkg; Pkg.develop(path=".");'
```

then run the following in the command line:
```
julia --project=docs -e 'using AlgebraicABMs, LiveServer; servedocs(literate_dir="docs/literate",skip_dir="docs/src/generated")'
```
generating the document templates will take some time.

To locally run the test suite, first install the test framework and register this package:

```
julia --project=test -e 'using Pkg; Pkg.develop(path="."); Pkg.add("Aqua")'
```

and then run the following command:
```
julia --project=test test/runtests.jl
```

## NOTE
This library is currently under active development, and so is not yet at a
point where a constant API/behavior can be assumed. That being said, if this
project looks interesting/relevant please contact us and
[let us know](https://www.algebraicjulia.org/#contributing)!

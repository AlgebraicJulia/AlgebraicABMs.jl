# AlgebraicABMs.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://AlgebraicJulia.github.io/AlgebraicABMs.jl/stable)
[![Development Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://AlgebraicJulia.github.io/AlgebraicABMs.jl/dev)
[![Code Coverage](https://codecov.io/gh/AlgebraicJulia/AlgebraicABMs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AlgebraicJulia/AlgebraicABMs.jl)
[![CI/CD](https://github.com/AlgebraicJulia/AlgebraicABMs.jl/actions/workflows/julia_ci.yml/badge.svg)](https://github.com/AlgebraicJulia/AlgebraicABMs.jl/actions/workflows/julia_ci.yml)

Important example files: `docs/literate/petri_example.jl` and `test/ABMs.jl`

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

### To-do: Buildkite

AlgebraicJulia uses [Buildkite](https://buildkite.com/) to submit resource-intensive processes such as building documentation and executing tests to the [HiPerGator](https://www.rc.ufl.edu/about/hipergator/) computing cluster.

While this template comes with a preconfigured `.buildkite/pipeline.yml` file, this repository is not integrated with Buildkite by default. If you would like your repository to use Buildkite to run processes on HiPerGator, tag an issue with @AlgebraicJulia/SysAdmins. 

### 📔 To-do: Set Up GitHub Pages (Public Repos Only)

1. Follow the Usage steps above to set up a new template, make sure all initial GitHub Actions have passed
2. Navigate to the repository settings and go to "Code and automation", "Pages"
3. Make sure the "Source" dropdown is set to "Deploy from a branch"
4. Set the "Branch" dropdown to "gh-pages", make sure the folder is set to "/ (root)", and click "Save"
5. Go back to the main page of your repository and click the gear to the right of the "About" section in the right side column
6. Under "Website" check the checkbox that says "Use your GitHub Pages website" and click "Save changes"
7. You will now see a URL in the "About" section that will link to your package's documentation

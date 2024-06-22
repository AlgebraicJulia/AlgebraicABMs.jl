using Documenter
using Literate

const LITERATE_INPUT = joinpath(@__DIR__, "literate")
const LITERATE_OUTPUT = joinpath(@__DIR__, "src", "generated")

@info "Loading AlgebraicABMs"
using AlgebraicABMs

const no_literate = "--no-literate" in ARGS
if !no_literate
  @info "Building Literate.jl docs"

  # Set Literate.jl config if not being compiled on recognized service.
  config = Dict{String,String}()
  if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
    config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/AlgebraicABMs.jl/blob/gh-pages/dev"
    config["repo_root_url"] = "https://github.com/AlgebraicJulia/AlgebraicABMs.jl/blob/main/docs"
  end

  # for (root, dirs, files) in walkdir(literate_dir)
  #   out_dir = joinpath(generated_dir, relpath(root, literate_dir))
  #   for file in files
  #     f, l = splitext(file)
  #     if l == ".jl" && !startswith(f, "_")
  #       Literate.markdown(joinpath(root, file), out_dir;
  #         config=config, documenter=true, credit=false)
  #       Literate.notebook(joinpath(root, file), out_dir;
  #         execute=true, documenter=true, credit=false)
  #     end
  #   end
  # end

  for (root, _, files) ∈ walkdir(LITERATE_INPUT), file ∈ files
    # ignore non julia files
    splitext(file)[2] == ".jl" || continue
    # full path to a literate script
    ipath = joinpath(root, file)
    # generated output path
    opath = splitdir(replace(ipath, LITERATE_INPUT=>LITERATE_OUTPUT))[1]
    # generate the markdown file calling Literate
    Literate.markdown(ipath, opath)
  end

end

@info "Building Documenter.jl docs"
makedocs(
  modules=[AlgebraicABMs],
  format=Documenter.HTML(size_threshold=typemax(Int)),
  sitename="AlgebraicABMs.jl",
  doctest=false,
  checkdocs=:none,
  pages=Any[
    "AlgebraicABMs.jl"=>"index.md",
    "Examples"=>Any[
      "generated/sir_petri.md",
      "generated/game_of_life.md",
      "generated/lotka_volterra.md",
    ],
    "Library Reference"=>"api.md",
  ]
)

@info "Deploying docs"
deploydocs(
  target="build",
  repo="github.com/AlgebraicJulia/AlgebraicABMs.jl.git",
  branch="gh-pages"
)

using BenchmarkTools

const SUITE = BenchmarkGroup()

module BenchmarkABMs
  include("ABMs.jl")
end

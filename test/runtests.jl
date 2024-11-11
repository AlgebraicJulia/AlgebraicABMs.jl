using Test

@testset "Code Quality (Aqua.jl)" begin
  include("aqua.jl")
end

@testset "ABMs" begin
  include("ABMs.jl")
end

@testset "PetriInterface" begin
  include("PetriInterface.jl")
end

@testset "StateChartsInterface" begin
  include("TestStateChartsInterface")
end

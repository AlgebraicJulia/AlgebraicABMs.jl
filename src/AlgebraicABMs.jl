""" Some description of ths package
"""
module AlgebraicABMs

export PetriNetCSet

using Reexport

include("Upstream.jl")
include("ABMs.jl")
include("Distributions.jl")
include("Visualization.jl")

@reexport using .Upstream
@reexport using .Distributions
@reexport using .ABMs
@reexport using .Visualization

# Methods to be implemented by extensions
function PetriNetCSet end 

end # module

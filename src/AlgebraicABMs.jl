""" Some description of ths package
"""
module AlgebraicABMs

export PetriNetCSet, StateChartABMSchema_MultipleObjects, make_infectious_rule_MultipleObjects, radomlyAssignInitialInfectives, make_ABM

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

# statecharts
function StateChartABMSchema_MultipleObjects end
function make_infectious_rule_MultipleObjects end
function radomlyAssignInitialInfectives end
function make_ABM end



end # module

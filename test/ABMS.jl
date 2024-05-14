module TestABMs

# ENV["JULIA_DEBUG"] = "AlgebraicABMs"

using Test
using AlgebraicABMs

using Catlab, AlgebraicRewriting

using AlgebraicABMs.ABMs: ABMRule, DiscreteHazard, ContinuousHazard, RegularP, 
                          EmptyP, RepresentableP, RuntimeABM

create_vertex = ABMRule(Rule(id(Graph()), create(ob(terminal(Graph)))), DiscreteHazard(1.))
@test create_vertex.pattern_type == EmptyP()

add_loop = ABMRule(Rule(id(Graph(1)), delete(Graph(1))), DiscreteHazard(1.5))
@test add_loop.pattern_type isa RepresentableP

rem_loop = ABMRule(Rule(delete(Graph(1)), id(Graph(1))), DiscreteHazard(2))
@test rem_loop.pattern_type == RegularP()

rem_edge = ABMRule(Rule(homomorphism(Graph(2), path_graph(Graph, 2); 
                        monic=true), id(Graph(2))), ContinuousHazard(1))
@test rem_edge.pattern_type isa RepresentableP


G = @acset Graph begin V=3; E=3; src=[1,1,1]; tgt=[1,1,2] end
abm = ABM([create_vertex, add_loop, rem_loop, rem_edge])

@test only(RuntimeABM(abm, G).clocks[3].val.match_vect)[1][:E](1) == 1 # One cached hom


traj = run!(abm, G);

end # module

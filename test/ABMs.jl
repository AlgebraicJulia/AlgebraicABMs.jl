module TestABMs

ENV["JULIA_DEBUG"] = "AlgebraicABMs" # turn on @debug messages for this package

using Test
using AlgebraicABMs
using Catlab, AlgebraicRewriting

using AlgebraicABMs.ABMs: RegularP, EmptyP, RepresentableP, RuntimeABM
using AlgebraicRewriting.Incremental.IncrementalCC: match_vect

# L = ∅, I = ∅, R = •↺
create_loop = ABMRule(
  :CreateLoop,
  Rule(id(Graph()),                  # l : I -> L
       create(ob(terminal(Graph)))), # r : I → R
  DiscreteHazard(1.)) # Dirac delta, indep. of clock time / state

# check that we know this rule has an empty pattern L
@test create_loop.pattern_type == EmptyP()

# • ← • → •↺
add_loop = ABMRule(Rule(id(Graph(1)), # 
                        delete(Graph(1))),  # r : I -> R
                   ContinuousHazard(1.5))
@test add_loop.pattern_type == RepresentableP(Dict(:V=>[1]))

# •↺ ⇽ • → •
rem_loop = ABMRule(:RemLoop, Rule(delete(Graph(1)), id(Graph(1))), DiscreteHazard(2))
@test rem_loop.pattern_type == RegularP()

# •→• ⇽ • → •
rem_edge = ABMRule(:RemEdge, Rule(homomorphism(Graph(2), path_graph(Graph, 2); initial=(V=1:2,)), 
                        id(Graph(2))), 
                   ContinuousHazard(1))
@test rem_edge.pattern_type == RepresentableP(Dict(:E=>[1]))

# Create initial state 
G = @acset Graph begin V=3; E=3; src=[1,1,1]; tgt=[1,1,2] end
to_graphviz(G)

# Assemble rules into ABM
abm = ABM([create_loop, add_loop, rem_loop, rem_edge])

# 2 loops, so 2 cached homs for the only rule with an explicit hom set
@test length(only(match_vect(RuntimeABM(abm, G)[:RemLoop].val))) == 2

traj = run!(abm, G; maxevent=10);

traj = run!(ABM([rem_edge]), G);
@test length(traj) == 3


traj = run!(ABM([add_loop]), G);
@test length(traj) > 3 # after we add a loop, the match persists + is resampled



# Test events in parallel
create_loop = ABMRule(Rule(id(Graph(1)), delete(Graph(1))), DiscreteHazard(1.));
create_vert = ABMRule(Rule(id(Graph()),  create(Graph(1))), DiscreteHazard(1.));
abm = ABM([create_loop, create_vert]);
traj = run!(abm, Graph(); maxtime=5);

# ODEs
######
using AlgebraicABMs, AlgebraicRewriting, Catlab

# State of world: a set of free-floating Float64s 
@present SchLSet(FreeSchema) begin X::Ob; D::AttrType; f::Attr(X, D) end 
@acset_type LSet(SchLSet){Float64} 

# Rule: copy a variable
v = @acset LSet begin X=1; D=1; f=[AttrVar(1)] end
v2 = @acset LSet begin X=2; D=1; f=[AttrVar(1), AttrVar(1)] end
dup_vertex = ABMRule(Rule(id(v), homomorphism(v, v2; initial=(X=[1],))), DiscreteHazard(1.))

# Dynamics: for an individual variable, it grows linearly w/ time
flow = ABMFlow(v, RawODE([_ -> 1.0]), :Grow, [], [(:D => 1)])
# Make ABM
abm = ABM([dup_vertex], [flow])

# Initial state
init = @acset LSet begin X=2; f=[1.1, 2.2] end

# res = run!(abm, init, maxevent=2)





end # module

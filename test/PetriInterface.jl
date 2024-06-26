module TestPetriInterface 

using Test
using AlgebraicABMs, Catlab
using AlgebraicABMs.ABMs: RepresentableP, RegularP
n = length(methods(ABM))
using AlgebraicPetri # check that package dependencies work
@test length(methods(ABM)) == n+1

using Distributions

ENV["JULIA_DEBUG"] = "" # turn off @debug messages

# ## Petri-net based model
#
# We define an SIRS model with birth and death

sir_pn= @acset LabelledPetriNet begin
  S=3; sname=[:S,:I,:R]
  T=7; tname=[:inf,:rec,:birth,:deathS,:deathI,:deathR,:wane]
  I=7; it=[1,1,2,4,5,6,7]; is=[1,2,2,1,2,3,3]
  O=5; ot=[1,1,2,3,7]; os=[2,2,3,1,1]
end;

# We declare parameters to specify the random waiting times
pS, pI, pR = 95,5,0
pop = 100
lifespan = 65*365
μ = 1/lifespan
β = 0.001
wane = 60;

# Create the ABM by associating stochastic timers with each transition
abm = ABM(sir_pn, (inf=ContinuousHazard(1 / β),
                   rec=ContinuousHazard(Weibull(weibullpar(30, 5)...)),
                   birth=ContinuousHazard(1 / (μ*pop)),
                   deathS=ContinuousHazard(1 / μ),
                   deathI=ContinuousHazard(1 / μ),
                   deathR=ContinuousHazard(1 / μ),
                   wane=ContinuousHazard(wane)));

@test abm.rules[1].pattern_type isa RepresentableP
@test abm.rules[2].pattern_type == RegularP()
@test nparts.(Ref(codom(Catlab.left(abm.rules[1].rule))), [:S,:I,:R]) == [1,1,0] # S, I

# Initial state
init = PetriNetCSet(sir_pn; S=pS, I=pI)

# Run the model
res = run!(abm, init; maxtime=1000);

end # module

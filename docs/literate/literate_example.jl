# # Code Example
#
# First we want to load our package with `using`

using AlgebraicABMs
using Catlab, AlgebraicPetri
using Distributions
using Makie, CairoMakie

using AlgebraicABMs.PetriInterface: PetriNetCSet, make_rules

# ## Petri-net based model
#
# We define the SIR model


sir_pn= @acset LabelledPetriNet begin
  S=3; sname=[:S,:I,:R]
  T=7; tname=[:inf,:rec,:birth,:deathS,:deathI,:deathR,:wane]
  I=7; it=[1,1,2,4,5,6,7]; is=[1,2,2,1,2,3,3]
  O=5; ot=[1,1,2,3,7]; os=[2,2,3,1,1]
end
t1, t2, t3, t4, t5, t6, t7 = make_rules(sir_pn)

@assert PetriNetCSet(sir_pn, [1,2,3]) == PetriNetCSet(sir_pn; I=2,R=3,S=1)


## Parameters to specify the random waiting times
pS, pI, pR = 95,5,0
pop = 100
lifespan = 65*365
μ = 1/lifespan
β = 0.001
wane = 60

## functions which take a time point and return a distribution of waiting times
clockdists = Dict{Symbol,Function}()

# the Exponential clocks (Markov)
clockdists[:inf] = (t) -> Exponential(1 / β)
clockdists[:birth] = (t) -> Exponential(1 / (μ*pop))
clockdists[:deathS] = (t) -> Exponential(1 / μ)
clockdists[:deathI] = (t) -> Exponential(1 / μ)
clockdists[:deathR] = (t) -> Exponential(1 / μ)
clockdists[:wane] = (t) -> Exponential(wane)

# the Weibull clock (non-Markov)
α, θ = weibullpar(30, 5)
clockdists[:rec] = (t) -> Weibull(α, θ)


# ------------------------------------------------------------------------------

count(acs::ACSet) = 
  NamedTuple(Dict([o => nparts(acs, o) for o in ob(acset_schema(acs))]))

init = PetriNetCSet(sir_pn, [100,5,0])

# simulate a few models
sirout = run!(sir_pn, clockdists, init; save=count, maxevent=2000)
X = first.(sirout)
SIR = [getindex.(last.(sirout), x) for x in 1:3]
f = Figure();
ax = Axis(f[1,1])
ln1, ln2, ln3 = lines!.(Ref(ax), Ref(X), SIR)
Legend(f[1, 2], [ln1,ln2,ln3], ["S", "I","R"])
f
Makie.save("figures/SIRSDiscretetrajectory.png", f, px_per_unit=1, size=(800,600))

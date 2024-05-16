# # Petri Net rewriting
#
# First we want to load our package with `using`

using AlgebraicABMs
using Catlab, AlgebraicPetri # for declaring model building blocks 
using Distributions # for defining hazard rates
using Makie, CairoMakie # visualization

# ## Petri-net based model
#
# We define an SIRS model with birth and death

sir_pn= @acset LabelledPetriNet begin
  S=3; sname=[:S,:I,:R]
  T=7; tname=[:inf,:rec,:birth,:deathS,:deathI,:deathR,:wane]
  I=7; it=[1,1,2,4,5,6,7]; is=[1,2,2,1,2,3,3]
  O=5; ot=[1,1,2,3,7]; os=[2,2,3,1,1]
end

to_graphviz(sir_pn)

# A "state of the world" in this model is just a finite set of Susceptible, Infected, and Recovered people. Thus we can specify a model with 3 integers.

init = PetriNetCSet(sir_pn, S=100, I=5); # Initial state


# We declare parameters to specify the random waiting times
pS, pI, pR = 95,5,0
pop = 100
lifespan = 65*365
μ = 1/lifespan
β = 0.001
wane = 60;

# We associate with each transition a function which takes a time point and return a distribution of waiting times
clockdists = Dict{Symbol,Function}();

## the Exponential clocks (Markov)
clockdists[:inf] = (t) -> Exponential(1 / β)
clockdists[:birth] = (t) -> Exponential(1 / (μ*pop))
clockdists[:deathS] = (t) -> Exponential(1 / μ)
clockdists[:deathI] = (t) -> Exponential(1 / μ)
clockdists[:deathR] = (t) -> Exponential(1 / μ)
clockdists[:wane] = (t) -> Exponential(wane)

## the Weibull clock (non-Markov)
α, θ = weibullpar(30, 5)
clockdists[:rec] = (t) -> Weibull(α, θ);

# We make a reporting function that extracts information from each time step.
count(acs::ACSet) = 
  NamedTuple(Dict([o => nparts(acs, o) for o in ob(acset_schema(acs))]));


# We simulate a model
sirout = run!(sir_pn, clockdists, init; save=count, maxevent=2000)
X = first.(sirout)
SIR = [getindex.(last.(sirout), x) for x in 1:3]
f = Figure()
Legend(f[1, 2], lines!.(Ref(Axis(f[1,1])), Ref(X), SIR), ["S", "I","R"])

## We can also save the figure to the filesystem.
isdir("figures") || mkdir("figures");
Makie.save("figures/SIRSDiscretetrajectory.png", f, px_per_unit=1, size=(800,600));

f


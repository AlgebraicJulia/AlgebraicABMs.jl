# # Petri-net based ABMs
# 
# First we want to load our packages with `using

using AlgebraicABMs, Catlab
using AlgebraicPetri
using Distributions, Makie, CairoMakie
ENV["JULIA_DEBUG"] = ""; # turn off @debug messages for this package

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
# Initial state
init = PetriNetCSet(sir_pn; S=pS, I=pI)

# Run the model
res = run!(abm, init; maxtime=2000);

# Plot results
Makie.plot(res; Dict(o=>X->nparts(X,o) for o in [:S,:I,:R])...)


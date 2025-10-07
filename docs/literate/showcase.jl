# # Showcasing AlgebraicABMs features
# 
# This file is an incomplete showcase of some less obvious features in 
# AlgebraicABMs. It is not an intro-tutorial (the other demos are more geared 
# towards that). To begin, we want to load our packages with `using

using AlgebraicABMs, Catlab, AlgebraicRewriting
ENV["JULIA_DEBUG"] = "AlgebraicABMs"; # turn off @debug messages for this package

# ## Basis 

# Rewrite rules have associated timers which are scheduled to fire in the future.
# For normal rules, at any point in time, there is exactly one such timer per 
# match morphism from the pattern $L$ of the rule into the present state of the 
# world, $X$. This is often what we want, though in some circumstances, we want the 
# 'per-what?' of the rule to be something different from the pattern. Thus we 
# have the ability to explicitly give a *basis* for the rule, given as a 
# map $B \rightarrowtail L$ into the rule's pattern. The rule will have timers 
# per every 'match' $B \rightarrow X$ (which is just a *partial* match 
# $L \nrightarrow X$). Once it's time to fire, the rest of the match is 
# completed at random (if possible, otherwise nothing happens).

# For example, we may want to associate the timers of an infection process with 
# simply the infected people, even though the pattern of an infection rule is 
# a suceptible + an infected person. E.g. "every infected person tries to infect
# *someone* once per day", which would not be affected by how many susceptible 
# people there are. Let's demonstrate this:

@present SchSI(FreeSchema) begin
  (S,I)::Ob
end

@acset_type SI(SchSI)

si   = @acset SI begin S=1; I=1 end 
i    = @acset SI begin I=1 end 
ii   = @acset SI begin I=2 end 
init = @acset SI begin S=10000;I=1 end 
inf_rule = Rule(homomorphism(i, si), homomorphism(i, ii; any=true))
basis_hom = homomorphism(i, si) # our actual basis is just a single infected

abm_rule = ABMRule(inf_rule, DiscreteHazard(1); basis = basis_hom)
abm = ABM([abm_rule])
run!(abm, init; maxevent=3);

# We can see that after the first timestep there is 1 infection, then 2 on the 
# next timestep, then four on the next one.

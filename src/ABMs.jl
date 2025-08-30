
module ABMs

export ABM, ABMRule, run!, DiscreteHazard, ContinuousHazard, FullClosure, 
       ClosureState, ClosureTime, RawODE, ABMFlow, filter, push!, copy, length

using Distributions, CompetingClocks, Random
using DataStructures: DefaultDict
using DifferentialEquations: ODEProblem
using StructEquality

using Catlab, AlgebraicRewriting
using AlgebraicRewriting.Incremental.Algorithms: connected_acset_components, pull_back
using AlgebraicRewriting.Rewrite.Migration: repr_dict
using Catlab.CategoricalAlgebra.Chase: extend_morphism_constraints
using AlgebraicRewriting.Rewrite.Utils: get_pmap, get_rmap, get_expr_binding_map
import Catlab: left, right
import AlgebraicRewriting: get_match, ruletype, addition!, deletion!, get_matches

import ..Upstream: pattern, pops!, IncHomSet_basis

# Timers
########
# Key abstractions for timers commonly include closures that can use or remember 
# 0, 1 or 2 of 1) morphism (ACSetTransformation) 2) time information to produce a 
# distribution for a Distribution from a hazard rate

"""
Something that can produce a ACSetTransformation × clocktime → hazard_rate
"""
# A timer is broadly the same as a clock
abstract type AbsTimer end

abstract type StateDependentTimer <: AbsTimer end

state_dep(t::AbsTimer) = t isa StateDependentTimer

"""
A closure which accepts a ACSetTransformation and returns a function of type
clocktime → hazard_rate
"""
struct FullClosure <: StateDependentTimer
  val::Function # ACSetTransformation → clocktime → hazard_rate
end

(c::FullClosure)(m::ACSetTransformation, t::Float64) = c.val(m,t)

"""
A closure which accepts a clocktime and returns a hazard_rate. This is a timer 
which cannot depend on the match data nor ACSet state.
"""
struct ClosureTime <: AbsTimer
  val::Function # clocktime → hazard_rate
end

(c::ClosureTime)(t::Float64) = c.val(t)

"""
A closure which accepts a match morphism and returns a hazard_rate. This is a 
timer which cannot depend on the absolute clock time.
"""
struct ClosureState <: StateDependentTimer
  val::Function # ACSetTransformation → hazard_rate
end

(c::ClosureState)(m::ACSetTransformation) = c.val(m)


# An Abstract Hazard is a stateless, time invariant type of timer that 
# simply has a distribution (contiuous or discrete) associated with it
abstract type AbsHazard <: AbsTimer end

@struct_hash_equal struct DiscreteHazard <: AbsHazard
  val::Distribution{Univariate, Discrete}
end

DiscreteHazard(t::Number) = DiscreteHazard(Dirac(t))

@struct_hash_equal struct ContinuousHazard <: AbsHazard
  val::Distribution{Univariate, Continuous}
end

"""Check if a hazard rate is a simple exponential"""
# We can optimize simple exponentials by reasoning about drawing from the first of n (at a rate of n times the individual rate) rather than separately drawing from each of the n
is_exp(h::ContinuousHazard) = h.val isa Distributions.Exponential
is_exp(h::AbsTimer) = false


# Constant hazard rate
ContinuousHazard(p::Number) = ContinuousHazard(Exponential(p))




# Rules 
#######
# ASK: Why not AbsPatternType?
abstract type PatternType end

# Trivial map -- automatically matches.
"""Empty patterns have (one) trivial pattern match"""
@struct_hash_equal struct EmptyP <: PatternType end

"""
Default case, where pattern matches should be found via (incremental) 
homomorphism search and represented explicitly, each with own events getting 
scheduled.
"""
#**** NB: I believe that this would motivate searching naively via hom search as
#**** Here, there is no well (this likely doesn't need to be done incrementally)

@struct_hash_equal struct RegularP <: PatternType end

# """
# Special case of homsearch where no backtracking is needed. The only nonempty
# sets in L are those for objects with no outgoing homs. There may be attributes,
# however, so at runtime we must filter the sets before picking random elements.
# E.g. for labeled set L = {:a, :a, AttrVar(1)} we randomly pick two elements with
# label :a and one arbitrary element.

# WARNING: this is only viable if the timer associated with the rewrite rule is
# symmteric with respect to the discrete parts.
# """
# @struct_hash_equal struct DiscreteP <: PatternType
#   parts::Dict{Symbol, Int}
# end

# Here, we use the Yoneda lemma, and the fact that we are randomly selecting a 
# match anyway Hom(Rep_Y,G) is isomorphic to G[Y], where Rep_Y is the 
# representable for Y.
"""
A pattern match from a coproduct of representables is just a choice of parts
in the codomain. E.g. matching L = •→• • •  is just a random choice of edge and
two random vertices.

The vector of ints refers to parts of L which are the counits of the left kan 
extensions that define the representables (usually this is just wherever the 
colimit leg sends 1, as there is often just one X part in the representable X).

WARNING: this is only viable if the timer associated with the rewrite rule is
symmteric with respect to the disjoint representables and has a simple
exponential timer.
"""
@struct_hash_equal struct RepresentableP <: PatternType
  parts::Dict{Symbol, Vector{Int}}
end

# Check if this is specific to incremental hom search
Base.keys(p::RepresentableP) = keys(p.parts)

# Not sure what this is for
multiplier(p::RepresentableP, X::ACSet) =
  prod(nparts(X, k)^length(v) for (k, v) in pairs(p.parts))

# ASK: Meaning of this not_monic
not_monic(b::Bool) = b === false 
not_monic(obs::AbstractVector{Symbol}) = isempty(obs)



"""
Analyze a pattern to find the most efficient pattern type for it.

Because ACSet types do not know their own equations, we may have to pass the 
schema as an argument in order to compute representables that would otherwise 
be infinite.

Even if the pattern is a coproduct of representables, we cannot use the 
efficient encoding unless the distribution is either an exponential 
(or a single dirac delta - not yet supported).
"""
# I believe that this applies just as much to the case of naive hom search
# NB:  Rule is an AlgebraicRewriting Rule
function pattern_type(r::Rule, is_exp::Bool)
  p = pattern(r)
  
  # Check empty case
  isempty(p) && return EmptyP()

  # Determine if pattern is a coproduct of representables
  if is_exp && isempty(r.conditions) && not_monic(r.monic)
    repr_loc = DefaultDict{Symbol, Vector{Int}}(() -> Int[])
    reprs = repr_dict(typeof(p))
    ccs, iso′ = connected_acset_components(p)
    iso = invert_iso(iso′)
    for cc_leg in legs(ccs)
      found = false
      for (o, (repr, i)) in pairs(reprs)
        α = isomorphism(repr, dom(cc_leg)) 
        if !isnothing(α)
          push!(repr_loc[o], iso[o](cc_leg[o](α[o](i))))
          found = true
          break
        end
      end
      found || break
    end
    length(ccs) == sum(length.(values(repr_loc))) && return RepresentableP(repr_loc)
  end

  # Determine if pattern is discrete
  # all(ob(S)) do o 
  #   nparts(p, o) == 0  || isempty(homs(S, from=o))
  # end && return DiscreteP(Dict(o => nparts(p, o) for o in ob(S)))

  return RegularP() # no special case found
end




# Hazard rates depend on pattern type

# get_hazard functions basically go through the job of turning various types of
# hazard (hazards that state dependent, time dependent, state & time dependent, and stateless-memoryless-time independent,... )
# into a univariant distribution.
# The following apply to all pattern types
# This applies the hazard to the appropriate information that it needs to do its job
get_hazard(::PatternType, m::ACSetTransformation, t::Float64, h::FullClosure) = h(m, t)

get_hazard(::PatternType, ::ACSetTransformation, t::Float64, h::ClosureTime) = h(t)

get_hazard(::PatternType, m::ACSetTransformation, ::Float64, h::ClosureState) = h(m)

# This just gets the value of the associated distribution
get_hazard(::PatternType, ::ACSetTransformation, ::Float64, h::AbsHazard) = h.val

# ASK: why is this limited to representables?
# I suspect that this is dealing with the fact that for stateless, memoryless, time-invariant processes, with representables
# we can just get the first of the hazards.  
# ASK: why the adjustment by "multiplier(r,X)"?
function get_hazard(r::RepresentableP, f::ACSetTransformation, ::Float64, 
                    h::ContinuousHazard) 
   err = "Representable patterns must have simple exponential rules"
   X = codom(f)
   # TODO: It may be that the below "/" should be "*"??  Please check
   is_exp(h) ? Exponential(h.val.θ/multiplier(r,X)) : error(err)
end

const Maybe{T} = Union{Nothing, T}






"""
A stochastic rewrite rule with a dependent hazard rate

A basis is a subobject of the pattern of the rule for which we want a timer 
per match. By default, the basis ↣ pattern map is just id(pattern).

NB: A "basis" corresponds to the map from "L_fix" in the Edinburgh discussions 
(See "Lfix limits what Changes Are of Concern in Pattern Matching" in 
AlgebraicABMs slides)

"""

@struct_hash_equal struct ABMRule
  rule::Rule
  timer::AbsTimer
  basis::Maybe{ACSetTransformation}
  name::Maybe{Symbol}
  pattern_type::PatternType
  ABMRule(r::Rule, t::AbsTimer; basis=nothing, name=nothing) = 
    new(r, t, basis, name, pattern_type(r, is_exp(t)))
end

# Give name as first arg rather than as kwarg
ABMRule(name::Maybe{Symbol}, r::Rule, t::AbsTimer; kw...) = 
  ABMRule(r, t; name, kw...)


# Get the AlgebraicRewriting rewriting rule associated with this ABMRule
getrule(r::ABMRule) = r.rule

Base.nameof(r::ABMRule) = r.name

pattern_type(r::ABMRule) = r.pattern_type

# The next 3 accessor functions basically get pieces of/information on the rewriting rule
pattern(r::ABMRule) = pattern(getrule(r))

# These get the "L" and "R" components of the rewriting span   L <- I -> R
left(r::ABMRule) = left(getrule(r))
right(r::ABMRule) = right(getrule(r))

ruletype(r::ABMRule) = ruletype(getrule(r))

basis(r::ABMRule) = r.basis

basis_pattern(r::ABMRule) = isnothing(r.basis) ? codom(left(r)) : dom(basis(r))

get_matches(r::ABMRule, args...; kw...) = 
  get_matches(getrule(r), args...; kw...)

# ASK: Is this to allow for lifting of the ABMRule via Data Migration Functors?
(F::Migrate)(r::ABMRule) = 
  ABMRule(F(r.rule), r.timer; basis=F(r.basis), name=r.name)


# I am assuming that this will remain unchanged
"""
A type which implements AbsDynamics must be able to compile to an ODE for some 
set of variables.
"""
abstract type AbsDynamics end 

"""Use raw Julia functions to define an ODE"""
@struct_hash_equal struct RawODE <: AbsDynamics 
  dynam::Vector{Function}
end

# ASK: Interpretation?
""" Continuous dynamics """
@struct_hash_equal struct ABMFlow 
  pat::ACSet
  dyn::AbsDynamics
  name::Maybe{Symbol}
  acs::Vector{Condition} # application conditions
  mapping::Vector{Pair{Symbol, Int}} # pair pat's variables w/ dyn quantities
end 






"""
An agent-based model.
"""*
# Key structure for an ABM -- rules, continuous dynamics, names (what are these?)
# ASK: Are the names associated with the variables in the patterns in the rewrite rule?
@struct_hash_equal struct ABM
  rules::Vector{ABMRule}
  dyn::Vector{ABMFlow}
  # A map from the name of a rule to its index in "rules".
  names::Dict{Symbol, Int}
  function ABM(rules, dyn=[]) 
    names = Dict(n=>i for (i,n) in enumerate(nameof.(rules)) if !isnothing(n))
    new(rules, dyn, names)
  end
end


#additions(abm::ABM) = right.(abm.rules)

# Migrate an ABM with an Data Migrations functor
(F::Migrate)(abm::ABM) = ABM(F.(abm.rules), abm.dyn)

# This is adding another dispatch option for "getindex" and "filter", both of
# which are defined in the Julia standard library ("Base")
Base.getindex(abm::ABM, i::Int) = abm.rules[i]
Base.getindex(abm::ABM, n::Symbol) = abm.rules[abm.names[n]]

# Build another ABM with rules satifying the predicate f
Base.filter(f, abm::ABM) = filter(f, abm.rules) |> ABM


# Add a rule to an ABM.
# ASK: Why is this needed?
function Base.push!(abm::ABM, r::ABMRule; overwrite=false)
  if haskey(abm.names, r.name)
    overwrite || error("The ABM already has a rule with this name, set overwrite=true to replace")
    abm.rules[abm.names[r.name]] = r
  else
    push!(abm.rules, r)
    # Associated with the name of this rule, record the index of this rule in the set of rules.
    abm.names[r.name] = length(abm.rules)
  end
  abm
end


# Shallow Duplication of ABMs
Base.copy(abm::ABM) = abm.rules |> copy |> ABM # shallow - rules have same pointers
# A notion of length of an ABM
Base.length(abm::ABM) = length(abm.rules)



# Here, we are handling the first to fire
# These are the 3 different ways one can refer to a hom (a particular match).  
#    In the case where pattern is a representable, no data is needed (we just randomly sample)
#    Given that I have a representation of the homset, might refer to the matches with a single number or a pair, but specific to incremental hom search
#  This is declaring the TYPE of the default_sampler -- this is one of the samplers 
#     supported in CompetingClocks (https://github.com/adolgert/CompetingClocks.jl/blob/main/src/sample/firsttofire.jl)

const default_sampler = FirstToFire{
  Union{Pair{Int, Nothing},   # non-explicit homset (Handling Empty pattern & representable pattern)
        Pair{Int, ACSetTransformation}},  # explicit multiple connected component homset
  Float64}

"""
Data structure for maintaining simulation information while running an ABM
"""
# ASK: Runtime data structures for an ABM
mutable struct RuntimeABM
  state::ACSet
  #   "clocks" is in 1-to-1 correspondence with rules, and tells us, for each rule, 
  #         how to get the next event for that rule.
  #   "clocks" IS still needed for naive hom search -- for each, we still need a way of sampling from this.

  tnow::Float64
  nevent::Int

  const sampler::SSA # stochastic simulation algorithm
  const rng::Distributions.AbstractRNG
  const names::Dict{Symbol, Int}
  const prob::ODEProblem
  const probmap::Vector{Pair{Symbol, Int}}
  const probdict::Dict{Symbol, Dict{Int, Int}}

  function RuntimeABM(abm::ABM, init::T; sampler=default_sampler) where T<:ACSet
    # Create the runtime
    names = Dict(r => i for (i, r) in enumerate(nameof.(abm.rules))
                 if !isnothing(r))

  # sampler() is creating a new sampler (an empty schedule)
  rt = new(init, 
             0., 0, sampler(), Random.RandomDevice(), names, 
             mk_prob(abm, init)...)

    # Initialize the firing queue
    for (iRule, patType) in enumerate(pattern_type.(abm.rules))
      # ASK: should we get rid of this?

      # The keys here are the way of referring to an event of this sort, and the 
      # value is the homomorphism.
      eventKeysMorphisms = if patType isa RegularP
        [h => h for h in homomorphisms(patType, state)]
      else
        if patType isa EmptyP || all(>(0), nparts.(Ref(init), keys(patType)))
          # Here, we have a homomorphism that is either for an empty pattern or a
          # representable pattern that has all of its matches present.
          # The key is "nothing" because we don't need to refer to it 
          [nothing => create(init)]
        else 
          []
        end
      end

      for (eventKey, morphism) in keysMorphisms
        # Get the hazard rate.  This takes care of upsampling.
        # haz is a JULIA DISTRIBUTION.
        hazDistr = get_hazard(patType, morphism, 0., abm.rules[iRule].timer)
        # ASK: Is this notion of a sampler still relevant for naive hom search?
        # This will do the sampling of this event! (i )
        # We will a key
        # "key" for us is the hom itself (ACSetTransformation)

        # THIS IS SETTING UP INITIAL EVENTS -- this sampler is sampling the "haz"
        # This enable! is a "competing clocks" function
        # "haz" is a distribution.  Competing Clocks will sample from that distribution
        # "enable!" could schedule an event at a particular time, but may wait until letter
        # conceptually , 'sampler" has a schedule of events associated with it, and this
        # schedules one.  Sampler has no concept of the rule -- it sees no difference between events.
        # The sampler associates the new event at  a particular time with with the data "iRule => key"
        # Each event has its own hazard rate (or, hazard distribution from which it samples to schedule it)
        # Sampler.  Has no notion of types of events.
        # When want to add something to queue, call enable!
        # When want to delete something from the queue and get it, call "pops!"
        # When want to delete something from the queue without getting it, call "disable!"
        # Thinik of sampler as a smart dictionary of map of key time-of-next-fire
        # The only thing that identifies what an event is "key".  We label "iRule" to 
        # prevent clashes between events for rules sharing the same key.
        enable!(rt.sampler, iRule => eventKey, hazDistr, 0., 0., rt.rng)
      end
    end
    return rt
  end
end

state(r::RuntimeABM) = r.state




# ASK: What is the function of this haskey mechanism?  Is this specific to incremental hom search?  Is this to check if things are scheduled?  Is k something that matches a Keytype?
# REMOVE any of the below?
Base.haskey(rt::RuntimeABM, k::Pair) = haskey(rt.sampler.transition_entry, k)

Base.haskey(rt::RuntimeABM, k::Int) = 
  haskey(rt.sampler.transition_entry, k => nothing)

#Base.getindex(rt::RuntimeABM, i::Int) = rt.clocks[i]
#Base.getindex(rt::RuntimeABM, n::Symbol) = rt.clocks[rt.names[n]]

"""
Construct an ODE for a given ACSet state. Return a mapping which allows to go from index to AttrType+index. 
"""
# ASK: To better understand the mix of continuous and discrete, I'd love to better understand this.
function mk_prob(abm::ABM, state::ACSet)
  isempty(abm.dyn) && return (ODEProblem((_,_,_,_)->0, 0, (0.,1.)), [], Dict())
  error("HERE")
end




# Naive Hom searches to look for patterns to see which apply AFTER a state update
# Presumably we have to do that AFTER the actual state update
"""Pop the next random event, advance the clock"""
function pops!(rt::RuntimeABM)::Vector{Pair{Int, Maybe{ACSetTransformation}}}
  rt.nevent += 1
   # "which" are all the keys for all of the events scheduled to fire at this new time.
   # pops! (defined in Upstream.jl) : """"Get all the next events which occur simultaneously and disable them"""
  (rt.tnow, which) = pops!(rt.sampler, rt.rng, rt.tnow)
  return which
end


# Note the reference to get_match, but passing a timer that is (oddly) a homSet
# note the explicit homomorphism search!
# ASK: Should we remove the timer here?  
# L and G refer to the corresponding quantities in DPO rewriting.

# Handle an the explicit basis, by delegating to the other get_match morphisms
function get_match(pat::PatternType, L::ACSet, G::ACSet, m::Maybe{ACSetTransformation}; 
                   basis::Maybe{ACSetTransformation})::ACSetTransformation
  if isnothing(basis)
      return get_match(pat, L, G, m)
  else
    # Handle an explicit basis
    m = get_match(pat, dom(basis), G, m)
    initial = extend_morphism_constraints(m, basis)
    
    # ASK: Why is this a random choice?
    rand(homomorphisms(L, G; initial))
  end
end


# TODO
"""
Get match returns a randomly chosen morphism for the aggregate rule
"""
# ASK: what is create(G) here?  Is this handling the search for representables, 
#    allowing each part to be drawn randomly from the appropriate Set for the 
#    "Head honcho" of the representable in the ACSet (per the Yoneda Lemma)?
# The below are seemingly handling get_match under different types of hom sets (nothing, representable, etc.)
get_match(::EmptyP, L::ACSet, G::ACSet, ::Nothing)::ACSetTransformation = create(G)

function get_match(P::RepresentableP, L::T, G::ACSet, 
                   ::Nothing)::ACSetTransformation where T<:ACSet
  # ASK: What is the logic behind this choice of initial?
  initial = Dict(map(collect(pairs(P.parts))) do (o, idxs) 
    o => Dict(idx => rand(parts(G, o)) for idx in idxs)
  end)
  # ASK: Does this return a random homomorphism if there are more than 1?
  return homomorphism(L, G; initial)
end

# Handle a regular pattern
get_match(::RegularP, L::ACSet, G::ACSet, ACSetTransformation m)::ACSetTransformation = m


"""
A trajectory of an ABM: each event time and result of `save`.
"""
# The structure that records the trajectory of the ABM
@struct_hash_equal struct Traj
  init::ACSet                                       # Presumably the initial state
  events::Vector{Tuple{Float64, Int, String, Any}}  # History of events -- String gives the rule name.  what are the Int and Any for?
  hist::Vector{Span{<:ACSet}}                       # Presumably this is a span of before-after?
end


Traj(x::ACSet) = Traj(x, Tuple{Float64, Int, String, Any}[], Span{ACSet}[])

# This gives a way to accumulate the trajectory.
function Base.push!(t::Traj, tup::Tuple{Float64,Int,String,Any,Span{<:ACSet}}) 
  (τ, rule, rulename, v, sp) = tup
  push!(t.events, (τ, rule, rulename, v))
  # ASK Help understand this: We're fine if the history is empty or if (I suspect) we are starting from the thing that we produced last time???
  isempty(t.hist) || codom(left(sp)) == codom(right(last(t.hist))) || error(
    "Bad history \n$(codom(left(sp))) \n!= \n$(codom(right(last(t.hist))))"
  )
  push!(t.hist, sp)
end


Base.isempty(t::Traj) = isempty(t.events)

Base.length(t::Traj) = length(t.events)

const MAXEVENT = 100


"""
Run an ABM, creating a fresh runtime + trajectory.

save - function applied to the ACSet state to produce the data that gets stored for every change in the model
dt - timestep for checking discrete events when running ODE dynamics.
"""
# This uses an ABM and creates a new RuntimeABM
# Note that "init"" is the initial state (T is a subtype of ACSet)
function run!(abm::ABM, init::T; save=_->nothing, maxevent=MAXEVENT, 
              maxtime=Inf, kw...) where T<:ACSet 
  run!(abm::ABM, RuntimeABM(abm, init; kw...), Traj(init); 
       save, maxtime, maxevent)
end

function run!(abm::ABM, rt::RuntimeABM, output::Traj;
              save=_->nothing, maxevent=MAXEVENT, maxtime=Inf, dt=0.1)
  maxevent = isinf(maxtime) ? maxevent : typemax(Int)
  # Helper functions that automatically incorporate the runtime `rt`
  getname(rule::Int)::String = 
    string(isnothing(abm.rules[rule].name) ? rule : abm.rules[rule].name)
  log!(rule::Int, sp::Span) = 
    push!(output, (rt.tnow, rule, getname(rule), save(rt.state), sp))
 
  # Main loop
  while rt.nevent < maxevent && rt.tnow < maxtime
    # TODO: isempty(abm.dyn) should be check that all flows sum to 0 
    if length(rt.sampler) == 0 && isempty(abm.dyn)
      @info "Stochastic scheduling algorithm ran out of events"
      return output
    end

    # ASK: Is this the next event that will trigger?  It seems that perhaps the first time is obtained from the sampler?
    new_time = first(next(rt.sampler, rt.tnow, rt.rng))
    if !isempty(abm.dyn) && dt < new_time 
        # ASK: presumably have to integrate forward the ODE here?  
      error("HERE")
    else
        # Get next event + unpack data -- all of these events are at a given time.
        # NB: These events can be associated with different rules

      # This is the events that occur simultaneously at the next timepoint.  
      events::Vector{Pair{Int,Maybe{ACSetTransformation}}} = pops!(rt) # updates the clock time
      
      # ASK: What is the significance of the length of the sampler?
      N = length(rt.sampler)


      # ASK: What is the significance of the length of the events?  The set of events that could go off at this time?
      # Determine if "Event" needs an s at its end due to the plural case
      s = length(events) > 1 ? "s" : ""
      # ASK: Are we prefering the first of the events in general, or just here as a convenience for printing?
      rname(e) = let r = first(e); n = abm.rules[r].name; isnothing(n) ? r : n end
      @debug ("Step $(length(output)): Event$s $(join(string.(rname.(events)), ", "))"
              *" | Fired @ t = $(round(rt.tnow, digits=2)) ($N queued)")

      # TODO some sort of check that the events are consistent with each other
      # or a randomization of their order

      update_data = [] # use to update incremental hom sets afterwards
      # execute all the events
      for (event, keySchedule) in events
        # rule' is the AlgebraicRewriting rule for "rule" !
        rule::ABMRule = abm.rules[event]
        rule′::Rule = getrule(rule)
        rule_type::Symbol = ruletype(rule)
        
        # If RegularPattern, we have an explicit match, otherwise randomly pick one
        m = get_match(pattern_type(rule), pattern(rule), rt.state, keySchedule; 
                      basis=basis(rule))

        # bring the match 'up to speed' given the previous (simultaneous) updates
        # ASK: Here we are taking the pullback of other matches for incremental hom sets only?  Or
        # ASK: Is this needed for e.g., game of life, for those that have already been performed?  Or that are all queuing up to fire?
        for (l, r) in first.(update_data)
          m = pull_back(l, m) ⋅ r
        end
        
        
        dpo = rule_type == :DPO ? (left(rule′), m) : nothing
        # check if dangling condition is satisfied
        # ASK: What is dangling condition?  This is where we have a rule which deletes a
        #       a vertex (and thus technically matches), but where we can't delete that because
        #       something deletes this. 
    
        # We will have to check if when ACTUALLY DELETING things (NOT invalidated earlier) -- 
        # needs to not have a link into this.  
        
        isnothing(dpo) || can_pushout_complement(ComposablePair(dpo...)) || continue
        # Excute rewrite rule and unpack results
        # See notes NDO slides on AlgebraicABMs

        rw_result = (rule_type, rewrite_match_maps(rule′, m))
        rmap_ = get_rmap(rw_result...)
        xmap = get_expr_binding_map(rule′, m, rw_result[2])
        (lft, rght_) = get_pmap(rw_result...)
        rmap, rght = compose.([rmap_,rght_], Ref(xmap))
        pmap = Span(lft, rght)
        rt.state = codom(rmap) # update runtime state
        log!(event, pmap)      # record event result

        # Remove this?, given the lack of need to support incremental hom sets?  But is this this needed somehow for simultaneous match rules more generally?
        # Update data will be all the events that fired.   
        push!(update_data, (pmap, rmap, dpo, right(rule′)))
      end
      
      # if no event at this time was actionable, due to dangling condition (see above).

      isempty(update_data) && continue 

      # TODO: This is some key work remaining
      # Ok, now go through and do the naive hom searches
      # for each ABMRule R
      #   Accumulate CURRENT_STATE_MATCHES of all the match morphisms mState for rule R found in the current state rt.state
      #   Loop through all events e that involve rule R in the sampler
      #       Check if event's ACSetHomormophisms or Representables or Empty e.mSchedule is in the IN CURRENT_STATE_MATCHES
      #         if e.mSchedule is not in CURRENT_STATE_MATCHES, delete! e from sampler (these are events that were created earlier due to matches to this rule, but no longer match anything in the state)
      #         if e.mSchedule is in the current state as name mState, we leave it there, and mark mState as having been found in CURRENT_STATE_MATCHES (these are events in the event queue that are still present in the state, and therefore still of interest)
      #   We go through all of CURRENT_STATE_MATCHES that are not marked as found in the scheduler, and add them to the schedule (these are the new matches for this rule)
      #if not, remove any events for this rule


      # If any of the matches that were fired are still preserved, re-enable
      # ASK: This logic might be suggesting incremental hom search reasoning -- how to modify for naive hom search?
      # TODO: adapt this to the naive hom search context
      for (event, key) in events
        if haskey(rt.clocks[event], key)
          enable!′(rt.clocks[event][key], event, key)
        end
      end
    end
  end
  return output
end

end # module
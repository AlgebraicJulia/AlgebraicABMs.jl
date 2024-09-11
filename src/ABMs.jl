
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

"""
Something that can produce a ACSetTransformation × clocktime → hazard_rate
"""
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

abstract type AbsHazard <: AbsTimer end

@struct_hash_equal struct DiscreteHazard <: AbsHazard
  val::Distribution{Univariate, Discrete}
end

DiscreteHazard(t::Number) = DiscreteHazard(Dirac(t))

@struct_hash_equal struct ContinuousHazard <: AbsHazard
  val::Distribution{Univariate, Continuous}
end

"""Check if a hazard rate is a simple exponential"""
is_exp(h::ContinuousHazard) = h.val isa Distributions.Exponential
is_exp(h::AbsTimer) = false

ContinuousHazard(p::Number) = ContinuousHazard(Exponential(p))


# Rules 
#######
abstract type PatternType end

"""Empty patterns have (one) trivial pattern match"""
@struct_hash_equal struct EmptyP <: PatternType end

"""
Default case, where pattern matches should be found via (incremental) 
homomorphism search and represented explicitly, each with own events getting 
scheduled.
"""
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

Base.keys(p::RepresentableP) = keys(p.parts)

multiplier(p::RepresentableP, X::ACSet) =
  prod(nparts(X, k)^length(v) for (k, v) in pairs(p.parts))

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

get_hazard(::PatternType, m::ACSetTransformation, t::Float64, h::FullClosure) = h(m, t)

get_hazard(::PatternType, ::ACSetTransformation, t::Float64, h::ClosureTime) = h(t)

get_hazard(::PatternType, m::ACSetTransformation, ::Float64, h::ClosureState) = h(m)

get_hazard(::PatternType, ::ACSetTransformation, ::Float64, h::AbsHazard) = h.val

function get_hazard(r::RepresentableP, f::ACSetTransformation, ::Float64, 
                    h::ContinuousHazard) 
   err = "Representable patterns must have simple exponential rules"
   X = codom(f)
   is_exp(h) ? Exponential(h.val.θ/multiplier(r,X)) : error(err)
end

const Maybe{T} = Union{Nothing, T}

"""
A stochastic rewrite rule with a dependent hazard rate

A basis is a subobject of the pattern of the rule for which we want a timer 
per match. By default, the basis ↣ pattern map is just id(pattern).

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

getrule(r::ABMRule) = r.rule

Base.nameof(r::ABMRule) = r.name

pattern_type(r::ABMRule) = r.pattern_type

pattern(r::ABMRule) = pattern(getrule(r))

left(r::ABMRule) = left(getrule(r))
right(r::ABMRule) = right(getrule(r))

ruletype(r::ABMRule) = ruletype(getrule(r))

basis(r::ABMRule) = r.basis

basis_pattern(r::ABMRule) = isnothing(r.basis) ? codom(left(r)) : dom(basis(r))

get_matches(r::ABMRule, args...; kw...) = 
  get_matches(getrule(r), args...; kw...)

(F::SimpleMigration)(r::ABMRule) = 
  ABMRule(F(r.rule), r.timer; basis=F(r.basis), name=r.name)

"""
A type which implements AbsDynamics must be able to compiled to an ODE for some 
set of variables.
"""
abstract type AbsDynamics end 

"""Use raw Julia functions to define an ODE"""
@struct_hash_equal struct RawODE <: AbsDynamics 
  dynam::Vector{Function}
end

""" Continuous dynamics """
@struct_hash_equal struct ABMFlow 
  pat::ACSet
  dyn::AbsDynamics
  name::Maybe{Symbol}
  acs::Vector{Condition} # application conditions
  mapping::Vector{Pair{Symbol, Int}} # pair pat's variables w/ dyn quantities
end 

# Accessing an IncHomSet
const KeyType = Union{Pair{Int, Int},        # connected comp. homset
                      Vector{Pair{Int,Int}}} # multi-component homset

"""
An agent-based model.
"""
@struct_hash_equal struct ABM
  rules::Vector{ABMRule}
  dyn::Vector{ABMFlow}
  names::Dict{Symbol, Int}
  function ABM(rules, dyn=[]) 
    names = Dict(n=>i for (i,n) in enumerate(nameof.(rules)) if !isnothing(n))
    new(rules, dyn, names)
  end
end

additions(abm::ABM) = right.(abm.rules)

(F::SimpleMigration)(abm::ABM) = ABM(F.(abm.rules), abm.dyn)

Base.getindex(abm::ABM, i::Int) = abm.rules[i]
Base.getindex(abm::ABM, n::Symbol) = abm.rules[abm.names[n]]

Base.filter(f, abm::ABM) = filter(f, abm.rules) |> ABM

function Base.push!(abm::ABM, r::ABMRule; overwrite=false)
  if haskey(abm.names, r.name)
    overwrite || error("The ABM already has a rule with this name, set overwrite=true to replace")
    abm.rules[abm.names[r.name]] = r
  else
    push!(abm.rules, r)
    abm.names[r.name] = length(abm.rules)
  end
  abm
end

Base.copy(abm::ABM) = abm.rules |> copy |> ABM # shallow - rules have same pointers
Base.length(abm::ABM) = length(abm.rules)

"""A collection of timers associated at runtime w/ an ABMRule"""
abstract type AbsHomSet end

@struct_hash_equal struct EmptyHomSet <: AbsHomSet end

@struct_hash_equal struct RepresentableHomSet <: AbsHomSet end

@struct_hash_equal struct ExplicitHomSet <: AbsHomSet val::IncHomSet end

Base.keys(h::ExplicitHomSet) = keys(h.val)

Base.haskey(h::ExplicitHomSet, k::KeyType) = haskey(h.val, k)

Base.haskey(::EmptyHomSet, k) = false

Base.haskey(::RepresentableHomSet, k) = false

Base.pairs(h::ExplicitHomSet) = pairs(h.val)

Base.getindex(h::ExplicitHomSet, i) = h.val[i]

deletion!(h::ExplicitHomSet, m; kw...) =  deletion!(h.val, m; kw...)

addition!(h::ExplicitHomSet, k, r, u) = addition!(h.val, k, r, u)

"""Initialize runtime hom-set given the rule and the initial state"""
function init_homset(rule::ABMRule, state::ACSet, 
                     additions::Vector{<:ACSetTransformation})
  p, sd = pattern_type(rule), state_dep(rule.timer)
  p == EmptyP() && return EmptyHomSet()
  (sd || p == RegularP()  
   ) && return ExplicitHomSet(IncHomSet_basis(getrule(rule), state,  additions; 
                                        basis=basis_pattern(rule)))
  @assert p isa RepresentableP  "$(typeof(p))"
  return RepresentableHomSet()
end 

const default_sampler = FirstToFire{
  Union{Pair{Int, Nothing},   # non-explicit homset
        Pair{Int, Pair{Int,Int}}, # explicit single cc homset
        Pair{Int, Vector{Pair{Int,Int}}}},  # explicit mc homset
  Float64}

"""
Data structure for maintaining simulation information while running an ABM
"""
mutable struct RuntimeABM
  state::ACSet
  const clocks::Vector{AbsHomSet}
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
    rt = new(init, init_homset.(abm.rules, Ref(init), Ref(additions(abm))), 
             0., 0, sampler(), Random.RandomDevice(), names, 
             mk_prob(abm, init)...)
    # Initialize the firing queue
    for (i, (pat,homset)) in enumerate(zip(pattern_type.(abm.rules), rt.clocks))
      kv = if homset isa ExplicitHomSet 
        pairs(homset) 
      else
        if pat isa EmptyP || all(>(0), nparts.(Ref(init), keys(pat)))
          [nothing => create(init)]
        else 
          []
        end
      end
      for (key, val) in kv
        haz = get_hazard(pat, val, 0., abm.rules[i].timer)
        enable!(rt.sampler, i => key, haz, 0., 0., rt.rng)
      end
    end
    return rt
  end
end

state(r::RuntimeABM) = r.state

Base.haskey(rt::RuntimeABM, k::Pair) = haskey(rt.sampler.transition_entry, k)

Base.haskey(rt::RuntimeABM, k::Int) = 
  haskey(rt.sampler.transition_entry, k => nothing)

Base.getindex(rt::RuntimeABM, i::Int) = rt.clocks[i]
Base.getindex(rt::RuntimeABM, n::Symbol) = rt.clocks[rt.names[n]]

"""
Construct an ODE for a given ACSet state. Return a mapping which allows to go from index to AttrType+index. 
"""
function mk_prob(abm::ABM, state::ACSet)
  isempty(abm.dyn) && return (ODEProblem((_,_,_,_)->0, 0, (0.,1.)), [], Dict())
  error("HERE")
end

"""
Check that RuntimeABM incremental hom sets have all valid homs.
"""
function validate(rt::RuntimeABM)
  for c in filter(c -> c isa IncHomSet, rt.clocks)
    c.state == rt.state || error("State mismatch")
    validate(c)
  end
end

"""Pop the next random event, advance the clock"""
function pops!(rt::RuntimeABM)::Vector{Pair{Int, Maybe{KeyType}}}
  rt.nevent += 1
  (rt.tnow, which) = pops!(rt.sampler, rt.rng, rt.tnow)
  return which
end


function get_match(pat::PatternType, L::ACSet, G::ACSet, timer::AbsHomSet, key; 
                   basis::Maybe{ACSetTransformation}) 
  isnothing(basis) && return get_match(pat, L, G, timer, key)
  m = get_match(pat, dom(basis), G, timer, key)
  initial = extend_morphism_constraints(m, basis)
  rand(homomorphisms(L, G; initial))
end

"""
Get match returns a randomly chosen morphism for the aggregate rule
"""
get_match(::EmptyP, L::ACSet, G::ACSet, ::EmptyHomSet, ::Nothing) = create(G)

function get_match(P::RepresentableP, L::T, G::ACSet, ::RepresentableHomSet, 
                   ::Nothing) where T<:ACSet
  initial = Dict(map(collect(pairs(P.parts))) do (o, idxs) 
    o => Dict(idx => rand(parts(G, o)) for idx in idxs)
  end)
  return homomorphism(L, G; initial)
end

get_match(::RegularP, ::ACSet, ::ACSet, hs::ExplicitHomSet, key::KeyType) = hs[key]

"""
A trajectory of an ABM: each event time and result of `save`.
"""
@struct_hash_equal struct Traj
  init::ACSet
  events::Vector{Tuple{Float64, Int, String, Any}}
  hist::Vector{Span{<:ACSet}}
end

Traj(x::ACSet) = Traj(x, Tuple{Float64, Int, String, Any}[], Span{ACSet}[])

function Base.push!(t::Traj, tup::Tuple{Float64,Int,String,Any,Span{<:ACSet}}) 
  (τ, rule, rulename, v, sp) = tup
  push!(t.events, (τ, rule, rulename, v))
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
  disable!′(key::Pair) = disable!(rt.sampler, key, rt.tnow)
  disable!′(i::Int) = disable!′(i => nothing)
  function enable!′(m::ACSetTransformation, rule_id::Int, key::Maybe{KeyType}=nothing) 
    rule = abm.rules[rule_id]
    haz = get_hazard(pattern_type(rule), m, rt.tnow, rule.timer)
    enable!(rt.sampler, rule_id => key, haz, rt.tnow, rt.tnow, rt.rng)
  end

  # Main loop
  while rt.nevent < maxevent && rt.tnow < maxtime
    # TODO: isempty(abm.dyn) should be check that all flows sum to 0 
    if length(rt.sampler) == 0 && isempty(abm.dyn)
      @info "Stochastic scheduling algorithm ran out of events"
      return output
    end

    new_time = first(next(rt.sampler, rt.tnow, rt.rng))
    if !isempty(abm.dyn) && dt < new_time 
      error("HERE")
    else
      # Get next event + unpack data
      events::Vector{Pair{Int,Maybe{KeyType}}} = pops!(rt) # updates the clock time
      N = length(rt.sampler)

      s = length(events) > 1 ? "s" : ""
      rname(e) = let r = first(e); n = abm.rules[r].name; isnothing(n) ? r : n end
      @debug ("Step $(length(output)): Event$s $(join(string.(rname.(events)), ", "))"
              *" | Fired @ t = $(round(rt.tnow, digits=2)) ($N queued)")

      # TODO some sort of check that the events are consistent with each other
      # or a randomization of their order

      update_data = [] # use to update incremental hom sets afterwards
      # execute all the events
      for (event, key) in events
        rule::ABMRule, clocks::AbsHomSet = abm.rules[event], rt.clocks[event]
        rule′::Rule, rule_type::Symbol = getrule(rule), ruletype(rule)
        # If RegularPattern, we have an explicit match, otherwise randomly pick one
        m = get_match(pattern_type(rule), pattern(rule), rt.state, clocks, key; 
                      basis=basis(rule))
        # bring the match 'up to speed' given the previous (simultanous) updates
        for (l, r) in first.(update_data)
          m = pull_back(l, m) ⋅ r
        end
        dpo = rule_type == :DPO ? (left(rule′), m) : nothing
        # check if dangling condition is satisfied
        isnothing(dpo) || can_pushout_complement(ComposablePair(dpo...)) || continue
        # Excute rewrite rule and unpack results
        rw_result = (rule_type, rewrite_match_maps(rule′, m))
        rmap_ = get_rmap(rw_result...)
        xmap = get_expr_binding_map(rule′, m, rw_result[2])
        (lft, rght_) = get_pmap(rw_result...)
        rmap, rght = compose.([rmap_,rght_], Ref(xmap))
        pmap = Span(lft, rght)
        rt.state = codom(rmap) # update runtime state
        log!(event, pmap)      # record event result
        push!(update_data, (pmap, rmap, dpo, right(rule′)))
      end
      
      # if no event at this time was actionable, due to dangling condition
      isempty(update_data) && continue 

      # All other rules can potentially update in response to the current event
      for (i, (ruleᵢ, clocksᵢ)) in enumerate(zip(abm.rules, rt.clocks))
        pt = pattern_type(ruleᵢ)
        if pt == EmptyP() && i ∈ first.(events)
          enable!′(create(rt.state), i)
        elseif pt == RegularP() # update explicit hom-set w/r/t span Xₙ ↩ • -> Xₙ₊₁
          for ((lft, rght), rmap, dpo, rule_right) in update_data
            del_invalid, del_new = deletion!(clocksᵢ, lft; dpo)

            for d in del_invalid # disable clocks which are invalidated
              (i=>d) ∈ events || disable!′(i => d) # (event,key) already diabled
            end

            for a in del_new
              enable!′(clocksᵢ[a], i, a) 
            end
            add_invalid, add_new = addition!(clocksᵢ, rule_right, rmap, rght)

            for d in add_invalid # disable clocks which are invalidated
              (i=>d) ∈ (events) || disable!′(i => d) # (event,key) already diabled
            end
            for a in add_new
              enable!′(clocksᵢ[a], i, a) 
            end
          end
        elseif pt isa RepresentableP
          relevant_obs = keys(pt)
          Xs = left(first(first(update_data))), right(first(last(update_data)))
          # we need to update current timer if # of parts has changed
          if i ∈ first.(events) && all(>(0), nparts.(Ref(rt.state), relevant_obs))
            enable!′(create(rt.state), i)
          elseif !all(ob -> allequal(nparts.(codom.(Xs), ob)), relevant_obs)
            currently_enabled = haskey(rt, i)
            currently_enabled && disable!′(i) # Disable if active
            # enable new timer if possible to apply rule
            if all(>(0), nparts.(Ref(rt.state), relevant_obs))
              enable!′(create(rt.state), i) 
            end
          end
        end
      end
      # If any of the matches that were fired are still preserved, re-enable
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
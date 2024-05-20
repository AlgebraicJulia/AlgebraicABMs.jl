
module ABMs

export ABM, ABMRule, run!, DiscreteHazard, ContinuousHazard, FullClosure, 
       ClosureState, ClosureTime

using Distributions, Fleck, Random
using DataStructures: DefaultDict
using StructEquality

using Catlab, AlgebraicRewriting
using AlgebraicPetri: AbstractReactionNet
using AlgebraicRewriting.Incremental.Algorithms: connected_acset_components
using AlgebraicRewriting.Rewrite.Utils: get_pmap, get_rmap, get_expr_binding_map
import Catlab: right
import AlgebraicRewriting: get_match, ruletype, addition!, deletion!, get_matches

using ..Upstream: Migrate′, repr_dict
import ..Upstream: pattern

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
  val::Function # clocktime → hazard_rate
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
is_exp(h::AbsHazard) = false

ContinuousHazard(p::Number) = ContinuousHazard(Exponential(p))

get_hazard(m::ACSetTransformation, t::Float64, h::FullClosure) = h(m, t)

get_hazard(::ACSetTransformation, t::Float64, h::ClosureTime) = h(t)

get_hazard(m::ACSetTransformation, ::Float64, h::ClosureState) = h(m)

get_hazard(::ACSetTransformation, ::Float64, h::AbsHazard) = h.val

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

"""
Analyze a pattern to find the most efficient pattern type for it.

Because ACSet types do not know their own equations, we may have to pass the 
schema as an argument in order to compute representables that would otherwise 
be infinite.

Even if the pattern is a coproduct of representables, we cannot use the 
efficient encoding unless the distribution is either an exponential 
(or a single dirac delta - not yet supported).
"""
function pattern_type(r::Rule, is_exp::Bool; schema=nothing)
  p = pattern(r)
  S = isnothing(schema) ? acset_schema(p) : schema
  
  # Check empty case
  isempty(p) && return EmptyP()

  # Determine if pattern is a coproduct of representables
  if is_exp
    repr_loc = DefaultDict{Symbol, Vector{Int}}(() -> Int[])
    reprs = repr_dict(typeof(p), S)
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

"""
A stochastic rewrite rule with a dependent hazard rate
"""
@struct_hash_equal struct ABMRule
  rule::Rule
  timer::AbsTimer 
  pattern_type::PatternType
  ABMRule(r::Rule, t::AbsTimer; schema=nothing) = 
    new(r, t, pattern_type(r, is_exp(t); schema))
end

getrule(r::ABMRule) = r.rule

pattern_type(r::ABMRule) = r.pattern_type

pattern(r::ABMRule) = pattern(getrule(r))

right(r::ABMRule) = right(getrule(r))

ruletype(r::ABMRule) = ruletype(getrule(r))

get_matches(r::ABMRule, args...; kw...) = 
  get_matches(getrule(r), args...; kw...)

(F::Migrate′)(r::ABMRule) = 
  ABMRule(F.F(r.rule), r.timer; schema=F.F.delta ? F.dom : F.codom)

"""
A type which implements AbsDynamics must be able to compiled to an ODE for some 
set of variables.
"""
abstract type AbsDynamics end 

"""Use a Petri Net with rates"""
@struct_hash_equal struct PetriDynamics <: AbsDynamics 
  val::AbstractReactionNet
end

"""Use a Stock Flow diagram (possibly this could be in a package extension)"""
@struct_hash_equal struct StockFlowDynamics <: AbsDynamics 
  val::Any # TODO: integrate with StockFlow.jl. 
end

"""Continuous dynamics"""
@struct_hash_equal struct ABMFlow 
  pat::ACSet
  dyn::AbsDynamics
  mapping::Vector{Pair{Symbol, Int}} # pair pat's variables w/ dyn quantities
end 

const Maybe{T} = Union{Nothing, T}

# Accessing an IncHomSet
const KeyType = Union{Pair{Int, Int}}       # connected comp. homset
                      Vector{Pair{Int,Int}} # multi-component homset

"""
An agent-based model.
"""
@struct_hash_equal struct ABM
  rules::Vector{ABMRule}
  dyn::Vector{ABMFlow}
  ABM(rules, dyn=[]) = new(rules, dyn)
end

additions(abm::ABM) = right.(abm.rules)

(F::Migrate′)(abm::ABM) = ABM(F.(abm.rules), abm.dyn)

"""A collection of timers associated at runtime w/ an ABMRule"""
abstract type AbsHomSet end

@struct_hash_equal struct EmptyHomSet <: AbsHomSet end

@struct_hash_equal struct DiscreteHomSet <: AbsHomSet end

@struct_hash_equal struct ExplicitHomSet <: AbsHomSet val::IncHomSet end

Base.keys(h::ExplicitHomSet) = keys(h.val)

Base.haskey(h::ExplicitHomSet, k::KeyType) = haskey(h.val, k)

Base.pairs(h::ExplicitHomSet) = pairs(h.val)

Base.getindex(h::ExplicitHomSet, i) = h.val[i]

deletion!(h::ExplicitHomSet, m; kw...) =  deletion!(h.val, m; kw...)

addition!(h::ExplicitHomSet, k, r, u) = addition!(h.val, k, r, u)

"""Initialize runtime hom-set given the rule and the initial state"""
function init_homset(rule::ABMRule, state::ACSet, additions::Vector{<:ACSetTransformation})
  p, sd = pattern_type(rule), state_dep(rule.timer)
  p == EmptyP() && return EmptyHomSet()
  (sd || p == RegularP()  
   ) && return ExplicitHomSet(IncHomSet(getrule(rule), state,  additions))
  @assert p isa RepresentableP  "$(typeof(p))"
  return DiscreteHomSet()
end 

const default_sampler = FirstToFire{
  Union{Pair{Int, Nothing},   # non-explicit homset
        Pair{Int, Pair{Int,Int}}, # explicit single cc homset
        Pair{Int, Vector{Pair{Int,Int}}}},  # explicit mc homset
  Float64}

"""
Data @struct_hash_equal structure for maintaining simulation information while running an ABM
"""
mutable struct RuntimeABM
  state::ACSet
  const clocks::Vector{AbsHomSet}
  tnow::Float64
  nevent::Int
  const sampler::SSA # stochastic simulation algorithm
  const rng::Distributions.AbstractRNG
  function RuntimeABM(abm::ABM, init::T; sampler=default_sampler) where T<:ACSet
    rt = new(init, init_homset.(abm.rules, Ref(init), Ref(additions(abm))), 
             0., 0, sampler(), Random.RandomDevice())
    for (i, homset) in enumerate(rt.clocks)
      kv = homset isa ExplicitHomSet ? pairs(homset) : [nothing => create(init)]
      for (key, val) in kv
        haz = get_hazard(val, 0., abm.rules[i].timer)
        enable!(rt.sampler, i => key, haz, 0., 0., rt.rng)
      end
    end
    rt
  end
end

state(r::RuntimeABM) = r.state

Base.haskey(rt::RuntimeABM, k::Pair) = haskey(rt.sampler.transition_entry, k)

Base.haskey(rt::RuntimeABM, k::Int) = 
  haskey(rt.sampler.transition_entry, k => nothing)

"""Pop the next random event, advance the clock"""
function Fleck.next(rt::RuntimeABM)
  rt.nevent += 1
  (rt.tnow, which) = pop!(rt.sampler, rt.rng, rt.tnow)
  return which
end


"""
Get match returns a randomly chosen morphism for the aggregate rule

TODO incorporate the number of possibilities as a multiplier for the rate
"""
get_match(::EmptyP, L::ACSet, G::ACSet, ::EmptyHomSet, ::Nothing) = create(G)

function get_match(P::RepresentableP, L::T, G::ACSet, ::DiscreteHomSet, 
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
  events::Vector{Tuple{Float64, Int, Any}}
  hist::Vector{Span{<:ACSet}}
end

Traj(x::ACSet) = Traj(x, Pair{Float64, Any}[], Span{ACSet}[])

function Base.push!(t::Traj, τ::Float64,rule::Int, v::Any, sp::Span{<:ACSet}) 
  push!(t.events, (τ, rule, v))
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
"""
function run!(abm::ABM, init::T; save=_->nothing, maxevent=MAXEVENT, 
              maxtime=Inf, kw...) where T<:ACSet 
  run!(abm::ABM, RuntimeABM(abm, init; kw...), Traj(init); save, maxtime, maxevent)
end

function run!(abm::ABM, rt::RuntimeABM, output::Traj;
              save=_->nothing, maxevent=MAXEVENT, maxtime=Inf)

  # Helper functions that automatically incorporate the runtime `rt`
  log!(rule::Int, sp::Span) = push!(output, rt.tnow, rule, save(rt.state), sp)
  disable!′(key::Pair) = disable!(rt.sampler, key, rt.tnow)
  disable!′(i::Int) = disable!′(i => nothing)
  function enable!′(m::ACSetTransformation, rule::Int, key=nothing) 
    haz = get_hazard(m, rt.tnow, abm.rules[rule].timer)
    enable!(rt.sampler, rule => key, haz, rt.tnow, rt.tnow, rt.rng)
  end

  # Main loop
  while rt.nevent < maxevent && rt.tnow < maxtime
    if length(rt.sampler) == 0 
      @info "Stochastic scheduling algorithm ran out of events"
      return output
    end

    # Get next event + unpack data
    event::Int, key::Maybe{KeyType} = next(rt) # updates the clock time
    rule::ABMRule, clocks::AbsHomSet = abm.rules[event], rt.clocks[event]
    rule′::Rule, rule_type::Symbol = getrule(rule), ruletype(rule)

    @debug "$(length(output)): event $event fired @ $(rt.tnow)"
    # If RegularPattern, we have an explicit match, otherwise randomly pick one
    m = get_match(pattern_type(rule), pattern(rule), rt.state, clocks, key)
    dpo = rule_type == :DPO ? (left(rule′), m) : nothing

    # Excute rewrite rule and unpack results
    rw_result = (rule_type, rewrite_match_maps(rule′, m))
    rmap_ = get_rmap(rw_result...)
    xmap = get_expr_binding_map(rule′, m, rw_result[2])
    (lft, rght_) = get_pmap(rw_result...)
    rmap, rght = compose.([rmap_,rght_], Ref(xmap))
    pmap = Span(lft, rght)
    rt.state = codom(rmap) # update runtime state
    log!(event, pmap)      # record event result

    # All other rules can potentially update in response to the current event
    for (i, (ruleᵢ, clocksᵢ)) in enumerate(zip(abm.rules, rt.clocks))
      pt = pattern_type(ruleᵢ)
      if pt == EmptyP()
        i == event && enable!′(create(rt.state), i)
      elseif pt == RegularP() # update explicit hom-set w/r/t span Xₙ ↩ • -> Xₙ₊₁

        del_invalid, del_new = deletion!(clocksᵢ, lft; dpo)
        for d in del_invalid # disable clocks which are invalidated
          (i,d)==(event,key) || disable!′(i => d) # (event,key) already diabled
        end
        for a in del_new; enable!′(clocksᵢ[a], i, a) end

        add_invalid, add_new = addition!(clocksᵢ, event, rmap, rght) # rght: R → Xₙ₊₁
        for d in add_invalid # disable clocks which are invalidated
          (i,d)==(event,key) || disable!′(i => d) # (event,key) already diabled
        end
        for a in add_new; enable!′(clocksᵢ[a], i, a) end

        # If the match that just fired is still preserved
        if i == event && haskey(clocksᵢ, key)
          enable!′(clocksᵢ[key], i, key)
        end

      elseif pt isa RepresentableP
        relevant_obs = keys(pt)
        # we need to update current timer if # of parts has changed
        if i == event && all(>(0), nparts.(Ref(rt.state), relevant_obs))
          enable!′(create(rt.state), i)
        elseif !all(ob -> allequal(nparts.(codom.(pmap), Ref(ob))), relevant_obs)
          currently_enabled = haskey(rt, i)
          currently_enabled && disable!′(i) # Disable if active
          # enable new timer if possible to apply rule
          if all(>(0), nparts.(Ref(rt.state), relevant_obs))
            enable!′(create(rt.state), i) 
          end
        end
      end
    end
  end
  return output
end


end # module
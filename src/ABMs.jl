
module ABMs

export ABM

using Distributions, Fleck, Random
using DataStructures: DefaultDict

using Catlab, AlgebraicRewriting
using AlgebraicPetri: AbstractReactionNet
using AlgebraicRewriting.Incremental: connected_acset_components, key_dict
using AlgebraicRewriting.Rewrite.Migration: pres_hash
import Catlab: acset_schema, right, is_isomorphic
import ..PetriInterface: run!
import AlgebraicRewriting: get_match
import AlgebraicRewriting.Incremental: addition!, deletion!

# Possibly upstream
###################
is_isomorphic(f::FinFunction) = is_monic(f) && is_epic(f)

pattern(r::Rule) = codom(left(r))

acset_schema(r::Rule) = acset_schema(pattern(r))

Base.pairs(h::IncHomSet) = [k => h[k] for k in keys(key_dict(h))]

"""
Extract data in representable cache as a dictionary.
This is an intermediate step in yoneda_cache that should be factored out.
"""
function repr_dict(T::Type, S=nothing; cache="cache")::Dict{Symbol, Tuple{ACSet, Int}}
  S = Presentation(isnothing(S) ? T : S)
  yoneda_cache(T, S; cache)
  cache_dir = joinpath(cache, "$(nameof(T))_$(pres_hash(S))")
  Dict(map(nameof.(generators(S, :Ob))) do name
    path, ipath = joinpath.(cache_dir, ["$name.json", "_id_$name.json"])
    name => (read_json_acset(T, path), parse(Int,open(io->read(io, String), ipath)))
  end)
end


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

struct DiscreteHazard <: AbsHazard
  val::Distribution{Univariate, Discrete}
end

DiscreteHazard(t::Number) = DiscreteHazard(Dirac(t))

struct ContinuousHazard <: AbsHazard
  val::Distribution{Univariate, Continuous}
end

ContinuousHazard(p::Number) = ContinuousHazard(Exponential(p))

get_hazard(m::ACSetTransformation, t::Float64, h::FullClosure) = h(m,t)
get_hazard(::ACSetTransformation, t::Float64, h::ClosureTime) = h(t)
get_hazard(m::ACSetTransformation, ::Float64, h::ClosureState) = h(m)
get_hazard(::ACSetTransformation, ::Float64, h::AbsHazard) = h.val

# Rules 
#######
abstract type PatternType end

"""Empty patterns have (one) trivial pattern match"""
struct EmptyP <: PatternType end

"""
Default case, where pattern matches should be found via (incremental) 
homomorphism search and represented explicitly, each with own events getting 
scheduled.
"""
struct RegularP <: PatternType end

# """
# Special case of homsearch where no backtracking is needed. The only nonempty
# sets in L are those for objects with no outgoing homs. There may be attributes,
# however, so at runtime we must filter the sets before picking random elements.
# E.g. for labeled set L = {:a, :a, AttrVar(1)} we randomly pick two elements with
# label :a and one arbitrary element.

# WARNING: this is only viable if the timer associated with the rewrite rule is
# symmteric with respect to the discrete parts.
# """
# struct DiscreteP <: PatternType
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
symmteric with respect to the disjoint representables.
"""
struct RepresentableP <: PatternType
  parts::Dict{Symbol, Vector{Int}}
end 
Base.keys(p::RepresentableP) = keys(p.parts)

"""
Analyze a pattern to find the most efficient pattern type for it.
"""
function pattern_type(r::Rule)
  p = pattern(r)
  S = acset_schema(p)
  
  # Check empty case
  isempty(p) && return EmptyP()

  # Determine if pattern is a coproduct of representables
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

  # Determine if pattern is discrete
  # all(ob(S)) do o 
  #   nparts(p, o) == 0  || isempty(homs(S, from=o))
  # end && return DiscreteP(Dict(o => nparts(p, o) for o in ob(S)))

  return RegularP() # no special case found
end

"""
A stochastic rewrite rule with a dependent hazard rate
"""
struct ABMRule
  rule::Rule
  timer::AbsTimer 
  pattern_type::PatternType
  ABMRule(r::Rule, t::AbsTimer) = new(r, t, pattern_type(r))
end
getrule(r::ABMRule) = r.rule
pattern_type(r::ABMRule) = r.pattern_type
pattern(r::ABMRule) = pattern(getrule(r))
right(r::ABMRule) = right(getrule(r))

abstract type AbsDynamics end 

"""Use a petri net with rates"""
struct PetriDynamics <: AbsDynamics 
  val::AbstractReactionNet
end

"""Continuous dynamics"""
struct ABMFlow 
  pat::ACSet
  dyn::AbsDynamics
  mapping::Vector{Pair{Symbol, Int}} # pair pat's variables  w/ dyn quantities
end 

"""
An agent-based model.
"""
struct ABM
  rules::Vector{ABMRule}
  dyn::Vector{ABMFlow}
  ABM(rules, dyn=[]) = new(rules, dyn)
end
additions(abm::ABM) = right.(abm.rules)

"""A collection of timers associated at runtime w/ an ABMRule"""
abstract type AbsHomSet end 
struct EmptyHomSet <: AbsHomSet end
struct DiscreteHomSet <: AbsHomSet end
struct ExplicitHomSet <: AbsHomSet val::IncHomSet end
Base.keys(h::ExplicitHomSet) = keys(h.val)
Base.pairs(h::ExplicitHomSet) = pairs(h.val)
Base.getindex(h::ExplicitHomSet, i) = h.val[i]
deletion!(h::ExplicitHomSet, m) =  deletion!(h.val, m)
addition!(h::ExplicitHomSet, k, r, u) = addition!(h.val, k, r, u)

"""Initialize runtime hom-set given the rule and the initial state"""
function init_homset(rule::ABMRule, state::ACSet, additions::Vector{<:ACSetTransformation})
  p, sd = pattern_type(getrule(rule)), state_dep(rule.timer)
  p == EmptyP() && return EmptyHomSet()
  (sd || p == RegularP()) && return ExplicitHomSet(IncHomSet(pattern(rule), additions, state))
  @assert p isa RepresentableP  "$(typeof(p))"
  return DiscreteHomSet()
end 

const KeyType = Union{Pair{Int,Int}} # connected component homset
                      Tuple{Int,Vector{Pair{Int,Int}}} # multi-component homset 
const default_sampler = FirstToFire{
  Union{Pair{Int, Nothing},   # non-explicit homset
        Pair{Int, KeyType}},  # explicit homset
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
Base.haskey(rt::RuntimeABM, k) = haskey(rt.sampler.transition_entry, k)

"""Pick the next random event, advance the clock"""
function Fleck.next(rt::RuntimeABM)
  rt.nevent += 1
  (rt.tnow, which) = next(rt.sampler, rt.tnow, rt.rng)
  which
end


"""
Get match returns a randomly chosen morphism for the aggregate rule

TODO incorporate the number of possibilities as a multiplier for the rate
"""
get_match(::EmptyP, L::ACSet, G::ACSet, ::EmptyHomSet, ::Nothing) = create(G)
function get_match(P::RepresentableP, L::T, G::ACSet, ::DiscreteHomSet, ::Nothing) where T<:ACSet
  initial = Dict(map(collect(pairs(P.parts))) do (o, idxs) 
    o => Dict(idx => rand(parts(G, o)) for idx in idxs)
  end)
  homomorphism(L, G; initial)
end
get_match(::RegularP, ::ACSet, ::ACSet, hs::ExplicitHomSet, key::KeyType) = hs[key]

"""
A trajectory of an ABM: each event time and result of `save`.
"""
struct Traj
  init::ACSet
  events::Vector{Tuple{Float64, Int, Any}}
  hist::Vector{Span{<:ACSet}}
end
Traj(x::ACSet) = Traj(x, Pair{Float64, Any}[], Span{ACSet}[])

function Base.push!(t::Traj, τ::Float64,rule::Int, v::Any, sp::Span{<:ACSet}) 
  push!(t.events, (τ, rule, v))
  isempty(t.hist) || codom(left(sp)) == codom(right(last(t.hist))) || error(
    "Bad history \n$(codom(left(sp))) \n!= \n$(right(last(t.hist)))"
  )
  push!(t.hist, sp)
end

Base.isempty(t::Traj) = isempty(t.events)

Base.length(t::Traj) = length(t.events)

const MAXEVENT = 10

"""
Run an ABM, creating a fresh trajectory.
"""
function run!(abm::ABM, init::T; save=deepcopy, maxevent=MAXEVENT, maxtime=Inf, 
              kw...) where T<:ACSet 
  run!(abm::ABM, RuntimeABM(abm, init; kw...); save, maxevent)
end

function run!(abm::ABM, rt::RuntimeABM, output::Union{Traj,Nothing}=nothing; 
              save=deepcopy, maxevent=MAXEVENT, maxtime=Inf)
  output = isnothing(output) ? Traj(rt.state) : output

  log!(rule::Int, sp::Span) = push!(output, rt.tnow, rule, save(rt.state), sp)

  disable!′(which) = disable!(rt.sampler, which, rt.tnow)

  function enable!′(m::ACSetTransformation, rule::Int, key=nothing) 
    haz = get_hazard(m, rt.tnow, abm.rules[rule].timer)
    enable!(rt.sampler, rule => key, haz, rt.tnow, rt.tnow, rt.rng)
  end

  while rt.nevent < maxevent && rt.tnow < maxtime
    which = next(rt) # get next event, update time
    isnothing(which) && return output # end b/c no more events!

    event, key = which
    rule, clock = abm.rules[event], rt.clocks[event]
    @debug "$(length(output)): event $event fired @ $(rt.tnow)"

    m = get_match(pattern_type(rule), pattern(rule), rt.state, clock, key)
    update_maps = rewrite_match_maps(getrule(abm.rules[event]), m)
    rh, kh, kg = update_maps[:rh], update_maps[:kh], update_maps[:kg]
    rt.state = codom(rh) # update runtime state
    log!(event, Span(kg, kh))  # record state after event

    if pattern_type(rule) == EmptyP()   # "always enabled" need special treatment
      disable!′(which)                  # their hom-set won't change, so clocks 
      enable!′(create(rt.state), event) # won't be reset, so do it manually
    end

    # update matches for all events
    for (t, (ruleₜ, clockₜ)) in enumerate(zip(abm.rules, rt.clocks))
      pt = pattern_type(ruleₜ)
      if pt == RegularP() # update explicit hom-set
        homs = clockₜ.val
        for d in Incremental.deletion!(clockₜ, kg)
          disable!′(t => d) # disable clocks which are invalidated
        end
        for a in Incremental.addition!(clockₜ, event, rh, kh)
          enable!′(clockₜ[a], t, a)
        end
      elseif pt isa RepresentableP
        if !all(o->all(is_isomorphic, [kg[o], kh[o]]), keys(pt))
          haskey(rt, t => nothing) && disable!′(t => nothing)
          if all(>(0), nparts.(Ref(rt.state), collect(keys(pt)))) 
            enable!′(create(rt.state), t)
          end
        end
      end
    end
  end
  return output
end


end # module
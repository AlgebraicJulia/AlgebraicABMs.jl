module RewriteSemiMarkov

export run!

using Catlab, AlgebraicRewriting, AlgebraicPetri
using Random
using Fleck
# --------------------------------------------------------------------------------
# we want something to store the rules, clocks associated to each, and their type

"""
  A presentation of a clock system, which is responsible for maintaining all the homsets for each
event possible in the system, and the set of enabled clocks for homset. It also stores the functions that
return a distribution over possible waiting times given a current simulation time.
"""
@present SchClockSystem(FreeSchema) begin
  Clock :: Ob # a single ID
  Event :: Ob # an entire class of clocks (e.g. "typed" clocks)
  event::Hom(Clock,Event)
  (NameType, RuleType, DistType, MatchType, KeyType, RngType, SampType)::AttrType 
  name::Attr(Event,NameType) # each event has a name
  rule::Attr(Event,RuleType) # each event has a rewrite rule
  dist::Attr(Event,DistType) # each event has a function that returns a firing time when it's enabled
  match::Attr(Event,MatchType) # each event has a homset
  key::Attr(Clock,KeyType)

  AlwaysEnabled :: Ob # a subset of Events
  always_enabled::Hom(AlwaysEnabled, Event)

  Global :: Ob # a set with one element, for ACSet-level parameters
  rng::Attr(Global, RngType) # random number generator
  sampler::Attr(Global, SampType) 
end


# const ClockKeyType = Tuple{Int,Int,Int} # for single connected component homsets
const ClockKeyType = Tuple{Int,Vector{Pair{Int,Int}}}

"""
The acset type that stores instances of `SchClockSystem`
"""
@acset_type AbsClockSystem(SchClockSystem, index=[:event], unique_index=[:name,:key])
const ClockSystem=AbsClockSystem{Symbol,Rule,Function,Incremental.IncSumHomSet,ClockKeyType, AbstractRNG, SSA}

rng(c::ClockSystem) = only(c[:rng])
sampl(c::ClockSystem) = only(c[:sampler])

"""
Given a `spn<:AbstractLabelledPetriNet` and a dictionary mapping each transition name to
a function that takes in parameter `t` and returns a distribution object, sample a single trajectory
of the stochastic dynamics on the petri net until `maxevent` events occur. Print stuff (or not) for debugging
with `verbose`.
"""
function run!(clocksys::ClockSystem, state::T; save=deepcopy, 
              maxevent=1000, verbose=false) where {T<:ACSet}
  nevent, tnow = 0, 0.0
  sample, RNG = sampl(clocksys), rng(clocksys)
  
  output = [0. => save(state)] # record initial state

  # when and what will happen next?
  (tnow, which) = next(sample, tnow, RNG)
  nevent += 1

  while nevent ≤ maxevent
    !verbose || println("event $(first(which)) fired at $tnow, total number of events: $(length(output))")
    event = first(which)
    update_maps = rewrite_match_maps(
      clocksys[event, :rule], 
      clocksys[event, :match][last(which)]
    )
    state = codom(update_maps[:rh])

    # record state after event
    push!(output, tnow => save(state))

    # "always enabled" transitions need special treatment when they fire
    # because their hom-set won't change, the clocks won't be reset, so we do it manually
    if !isempty(incident(clocksys, event, :always_enabled))
      disable!(sample, which, tnow)
      enable!(sample, which, clocksys[event, :dist](tnow), tnow, tnow, RNG)
    end    

    # update matches for all events
    for t in parts(clocksys, :Event)
      # update the hom-set
      del = Incremental.deletion!(clocksys[t,:match], update_maps[:kg])
      add = Incremental.addition!(clocksys[t,:match], event, update_maps[:rh], update_maps[:kh])
      del = map(del) do k
        (t,k)
      end
      add = map(add) do k
        (t,k)
      end
      # clocks that are disabled in the new marking (state):
      # 1: disable in the sampler
      # 2: remove from the clock system
      for c in del
        disable!(sample, c, tnow)
      end    
      rem_parts!(clocksys, :Clock, sort(vcat(incident(clocksys, del, :key)...)))
      # clocks that are newly enabled in the new marking (state):
      # 1: enable in the sampler
      # 2: add to the clock system
      for c in add
        enable!(sample, c, clocksys[t, :dist](tnow), tnow, tnow, RNG)
      end
      add_parts!(clocksys, :Clock, length(add), key = add, event = t)       
      
    end
    (tnow, which) = next(sample, tnow, RNG)
    nevent += 1
  end

  return output
end

end # module

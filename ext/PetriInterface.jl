module PetriInterface

using AlgebraicPetri
using Catlab, AlgebraicRewriting

using AlgebraicABMs
using AlgebraicABMs.ABMs: AbsDynamics
import AlgebraicABMs.ABMs: ABM
import AlgebraicABMs: PetriNetCSet

using StructEquality

function ABM(p::AbstractPetriNet, timers) 
  return ABM(map(parts(p, :T)) do t
    tname = transition_name(p, t)
    timer = haskey(timers, tname) ? timers[tname] : timers[t]
    ABMRule(tname, make_rule(p, t), timer)
  end, [])
end

"""Give a name for a Petri Net transition (name = label)"""
ob_name(pn::LabelledPetriNet, s::Int)::Symbol = pn[s, :sname]

"""Give a name for a Petri Net transition (e.g. name = "S3" for species #3)"""
ob_name(::PetriNet, s::Int)::Symbol = Symbol("S$s")

transition_name(pn::LabelledPetriNet, s::Int)::Symbol = pn[s, :tname]

transition_name(::PetriNet, s::Int)::Symbol = Symbol("T$s")


"""
Creates a discrete C-Set from a Petri net with one object for each species in
the Petri net. By default, this creates an *empty* C-Set instance, but there are
two ways one may also wish to specify how many tokens are in each species. One
can give a vector, where the indices correspond to the indices of the S table of
the petri net. Alternatively, one can give keyword arguments where the keys are
the names of the species (as determined by `ob_name`). 

For example, `PetriNetCSet(sir_labeled_pn, S=20, I=1)` would create an
*instance* of a C-Set that has three tables ("S","I","R"), no morphisms nor
attributes, and that instance would have 20 rows in the "S" table and 1 row in
the "I" table. In general, instances on this schema are effectively named tuples
(S::Int,I::Int,R::Int).
"""
function PetriNetCSet(pn::AbstractPetriNet, args=[]; kw...)
  res = AnonACSet(BasicSchema(ob_name.(Ref(pn), parts(pn, :S)),[]))
  for (arg, s) in zip(args, parts(pn, :S))
    add_parts!(res, ob_name(pn, s), arg)  # Add tokens from Vector{Int}
  end
  for (s, arg) in pairs(kw)
    add_parts!(res, s, arg)  # Add tokens by name
  end
  res
end

"""Assumes that tokens are deleted and recreated, rather than preserved"""
function make_rule(pn::Union{PetriNet, LabelledPetriNet}, t::Int)
  L, R = LR = [PetriNetCSet(pn) for _ in 1:2]
  add_part!.(Ref(L), ob_name.(Ref(pn), pn[incident(pn, t, :it), :is]))
  add_part!.(Ref(R), ob_name.(Ref(pn), pn[incident(pn, t, :ot), :os]))
  Rule(create.(LR)...)
end


"""Use a Petri Net with rates"""
@struct_hash_equal struct PetriDynamics <: AbsDynamics 
  val::AbstractReactionNet
end


end # module


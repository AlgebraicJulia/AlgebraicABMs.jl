"""For things that should likely be upstreamed, such as type piracy"""
module Upstream 

using Catlab, AlgebraicRewriting
import Catlab: is_isomorphic, Presentation
using AlgebraicRewriting.Rewrite.Migration: pres_hash
using Fleck: FirstToFire, disable!, next
using Distributions: AbstractRNG

# Upstream to Catlab
####################
Presentation(p::Presentation) = p

"""
Turn any span into a partial map by quotienting I and R by the left map: epi-mono factorize I → L and then take a pushout. 

                          I → R
                        ↙ ↡   ↓
                      L ↢ I' →⌜R'

Something like this may be needed if a Σ data migration fails to preserve monos for a DPO rewrite rule.
"""
function make_partial(s::Span{<:ACSet})
  i_i′, i′_l = epi_mono(left(s))
  i′_r, _  = pushout(i_i′, right(s))
  return Span(i′_l, i′_r)
end

# Upstream to AlgRewriting
##########################

# Fleck
#######
"""Get the next event and disable it. This will """
function Base.pop!(sampler::FirstToFire{K,T}, rng::AbstractRNG, 
                   tnow::T)::Tuple{T,K} where {K,T}
  isempty(sampler.firing_queue) && error("Cannot pop! an empty sampler")
  (new_time, which) = next(sampler, tnow, rng)
  disable!(sampler, which, new_time)
  return (new_time, which)
end 

"""Get all the next events which occur simultaneously and disable them"""
function pops!(sampler::FirstToFire{K,T}, rng::AbstractRNG, tnow::T
              )::Tuple{T,Vector{K}} where {K,T}
  (new_time, which) = pop!(sampler, rng, tnow)
  whiches = K[which]
  while true
    time, _ = next(sampler, tnow, rng)
    (time != new_time || isempty(sampler.firing_queue)) && break
    push!(whiches, last(pop!(sampler, rng, tnow)))
  end
  return (new_time, whiches)
end 


end # module

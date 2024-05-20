"""For things that should likely be upstreamed, such as type piracy"""
module Upstream 
export Migrate′

using Catlab, AlgebraicRewriting
import Catlab: is_isomorphic, Presentation
using AlgebraicRewriting.Rewrite.Migration: pres_hash
using Fleck: SSA, disable!, next
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

"""
TODO: extend AlgebraicRewriting.Migrate to include schema information about domain 
and codomain.
"""
struct Migrate′
  F::Migrate
  dom::Presentation
  codom::Presentation
  Migrate′(o::AbstractDict,
           h::AbstractDict,
           s1::Presentation,
           t1::Type,
           s2=nothing,
           t2=nothing; 
           delta::Bool=true) = 
    new(Migrate(o, h, t1, t2; delta), s1, isnothing(s2) ? s1 : s2)
end

(f::Migrate′)(x::ACSet) = f.F(x)
(f::Migrate′)(x::ACSetTransformation) = f.F(x)

# Fleck
#######
"""Get the next event and disable it"""
function Base.pop!(sampler::SSA{K,T}, rng::AbstractRNG, 
                   tnow::T)::Tuple{T,K} where {K,T}
  (new_time, which) = next(sampler, tnow, rng)
  disable!(sampler, which, new_time)
  return (new_time, which)
end 

end # module

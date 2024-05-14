# # Petri Net rewriting
#
# First we want to load our package with `using`

using AlgebraicABMs
using Catlab, AlgebraicPetri # for declaring model building blocks 
using Distributions # for defining hazard rates

# ## Schema 
# 
# We create a regular grid

@present SchLife(FreeSchema) begin 
  Cell::Ob 
  (N,E,W,S)::Hom(Cell,Cell)
  Life::Ob
  live::Hom(Life,Cell)
end

# We define a network of cells

@present SchLifeGraph <: SchSymmetricGraph begin 
  Life::Ob
  live::Hom(Life,V)
end

@acset_type Life(SchLifeGraph) <: AbstractSymmetricGraph

function living_neighbors(n::Int; alive=false)
  X = Life(1)
  alive && add_part!(X, :Live, live=1)
  for _ in 1:n
    v = add_part!(X, :V)
    add_part!(X, :Life, live=v)
    add_edge!(X, v, 1)
  end
  X
end

living_neighbors(2)


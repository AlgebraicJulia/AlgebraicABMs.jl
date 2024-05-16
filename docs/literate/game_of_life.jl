# # Game of Life
#
# First we want to load our package with `using`

using AlgebraicABMs, Catlab, AlgebraicRewriting

# ## Schema 
# 
# We define a network of cells that can be alive or dead (alive cells are in 
# the image of the `live` function, which picks out a subset of the vertices.)

@present SchLifeGraph <: SchSymmetricGraph begin 
  Life::Ob
  live::Hom(Life,V)
end

@acset_type Life(SchLifeGraph) <: AbstractSymmetricGraph;

to_graphviz(SchLifeGraph)

# We extend our schema with coordinate information - this doesn't affect the 
# game logic, so we  write our rules in the language of `Life` and only migrate 
# to `LifeCoords` at the end.

@present SchLifeCoords <: SchLifeGraph begin
  Coords::AttrType
  coords::Attr(V, Coords)
end

@acset_type LifeCoords(SchLifeCoords){Tuple{Int,Int}} <: AbstractSymmetricGraph;

to_graphviz(SchLifeCoords)

# Create a regular grid with periodic boundary conditions with an initial state
function make_grid(curr::AbstractMatrix)
  n, m = size(curr)
  n == m || error("Must be square")
  X, coords = LifeCoords(), Dict()
  for i in 1:n, j in 1:n
    coords[i=>j] = add_vertex!(X; coords=(i, j))
    Bool(curr[i, j]) && add_part!(X, :Life, live=coords[i=>j])
  end
  for i in 1:n, j in 1:n
    i < n && add_edge!(X, coords[i=>j], coords[i+1=>j])
    j < n && add_edge!(X, coords[i=>j], coords[i=>j+1])
    i < n && j < n && add_edge!(X, coords[i=>j], coords[i+1=>j+1])
    i < n && j > 1 && add_edge!(X, coords[i=>j], coords[i+1=>j-1])
  end
  return X
end;

# Create a random game state on a square, periodic grid
make_grid(n::Int, random=true) = 
  make_grid((random ? rand : zeros)(Bool, (n, n)));

# Visualize a game state (using coordinates if available)
function view_life(X::Union{Life, LifeCoords}, pth=tempname())
  pg = PropertyGraph{Any}(; prog="neato", graph=Dict(),
    node=Dict(:shape => "circle", :style => "filled", :margin => "0"),
    edge=Dict(:dir => "none", :minlen => "1"))
  add_vertices!(pg, nparts(X, :V))
  for v in vertices(X)
    set_vprop!(pg, v, :fillcolor, isempty(incident(X, v, :live)) ? "red" : "green")
    if X isa LifeCoords 
      x, y = X[v, :coords]
      set_vprop!(pg, v, :pos, "$x,$(y)!")
    end
  end
  for e in filter(e -> X[e, :inv] > e, edges(X))
    add_edge!(pg, X[e, :src], X[e, :tgt])
  end
  G = to_graphviz(pg)
  open(pth, "w") do io
    show(io, "image/svg+xml", G)
  end
  G
end;

view_life(make_grid(3)) # For example, a random 3x3 grid

    
#=
We ought to be able to take a state of the world (with no coordinate information)
and obtain a state of the world with coordinates (the canonical way to do this 
is to assign "variables" for the values of the coordinates).
=#

idₒ = Dict(x => x for x in Symbol.(generators(SchLifeGraph, :Ob)))
idₘ = Dict(x => x for x in Symbol.(generators(SchLifeGraph, :Hom)))
AddCoords = Migrate′(idₒ, idₘ, SchLifeGraph, Life, SchLifeCoords, LifeCoords; delta=false);
RemCoords = DeltaMigration(FinFunctor(idₒ, idₘ, SchLifeGraph, SchLifeCoords))
# ## Helper constants and functions 
const Dead = Life(1) # a single dead cell
const Live = @acset Life begin V=1; Life=1; live=1 end # a single living cell
const to_life = homomorphism(Dead, Live)  # the unique map Dead → Live

"""Create a context of n living neighbors for either a dead or alive cell"""
function living_neighbors(n::Int; alive=true)::ACSetTransformation
  X = Life(1)
  alive && add_part!(X, :Life, live=1)
  for _ in 1:n
    v = add_part!(X, :V)
    add_part!(X, :Life, live=v)
    add_edge!(X, v, 1)
  end
  homomorphism(alive ? Live : Dead, X; initial=(V=[1],))
end;

PAC(m) = AppCond(m; monic=true) # Positive Application condition
NAC(m) = AppCond(m, false; monic=true) # Negative Application condition
TickRule(args...; kw...) = # Rule which fires on 1.0, 2.0, ...
  ABMRule(Rule(args...; kw...), DiscreteHazard(1); schema=SchLifeGraph);

# ## Create model by defining update rules

# A cell is born iff it has three living neighbors
birth = TickRule(id(Dead), to_life; 
                 ac=[PAC(living_neighbors(3; alive=false)),
                     NAC(living_neighbors(4; alive=false)),
                     NAC(to_life)]);

# A cell is born iff it has ≥ 2 living neighbors but < 4 living neighbors
death = TickRule(to_life, id(Dead); 
                 ac=[PAC(living_neighbors(2)), 
                     NAC(living_neighbors(4))]);
 
GoL = ABM([birth, death]);


# Create an initial state
G = make_grid([1 0 1 0 1; 
               0 1 0 1 0; 
               0 1 0 1 0; 
               1 0 1 0 1;
               1 0 1 0 1])
view_life(G);
G = make_grid(ones(1,1))
# Run the model
migrate(Life, G, RemCoords)
get_matches(birth.rule, migrate(Life, G, RemCoords))
res = run!(AddCoords(GoL), G; maxevent=2);

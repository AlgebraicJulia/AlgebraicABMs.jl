# # Game of Life
#
# First we want to load our package with `using`

using AlgebraicABMs, Catlab, AlgebraicRewriting
ENV["JULIA_DEBUG"] = "AlgebraicABMs"; # turn on @debug messages for this package

# ## Schema 
# 
# We define a network of cells that can be alive or dead (alive cells are in 
# the image of the `live` function, which picks out a subset of the vertices.)

@present SchLifeGraph <: SchSymmetricGraph begin 
  Life::Ob
  live::Hom(Life,V)
end

# Create a datatype, 
@acset_type LifeState(SchLifeGraph) <: AbstractSymmetricGraph;

to_graphviz(SchLifeGraph) # visualize the schema

# We extend our schema with coordinate information - this doesn't affect the 
# game logic, so we  write our rules in the language of `Life` and only migrate 
# to `LifeCoords` at the end.

@present SchLifeCoords <: SchLifeGraph begin
  Coords::AttrType
  coords::Attr(V, Coords)
end

@acset_type LifeStateCoords(SchLifeCoords){Tuple{Int,Int}} <: AbstractSymmetricGraph;

to_graphviz(SchLifeCoords)

# Create a regular grid with empty boundary conditions with an initial state
function make_grid(curr::AbstractMatrix)
  n, m = size(curr)
  n == m || error("Must be square")
  X, coords = LifeStateCoords(), Dict()
  for j in 1:n, i in 1:n
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

# Create a random game state on a square grid
make_grid(n::Int, random=true) = 
  make_grid((random ? rand : zeros)(Bool, (n, n)));

# Visualize a game state (using coordinates if available)
function view_life(X::Union{LifeState, LifeStateCoords}, pth=tempname())
  pg = PropertyGraph{Any}(; prog="neato", graph=Dict(),
    node=Dict(:shape => "circle", :style => "filled", :margin => "0"),
    edge=Dict(:dir => "none", :minlen => "1"))
  add_vertices!(pg, nparts(X, :V))
  for v in vertices(X)
    set_vprop!(pg, v, :fillcolor, isempty(incident(X, v, :live)) ? "red" : "green")
    if X isa LifeStateCoords 
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
AddCoords = Migrate′(idₒ, idₘ, SchLifeGraph, LifeState, SchLifeCoords, LifeStateCoords; delta=false);
RemCoords = DeltaMigration(FinFunctor(idₒ, idₘ, SchLifeGraph, SchLifeCoords));
# ## Helper constants and functions 
const DeadCell = LifeState(1) # a single dead cell
const LiveCell = @acset LifeState begin V=1; Life=1; live=1 end # a single living cell
const to_life = homomorphism(DeadCell, LiveCell)  # the unique map Dead → Live

PAC(m) = AppCond(m; monic=true) # Positive Application condition
NAC(m) = AppCond(m, false; monic=true) # Negative Application condition
TickRule(name, args...; kw...) = # Rule which fires on 1.0, 2.0, ...
  ABMRule(name, Rule(args...; kw...), DiscreteHazard(1); schema=SchLifeGraph);

"""Create a context of n living neighbors for either a dead or alive cell"""
function living_neighbors(n::Int; alive=true)::ACSetTransformation
  X = LifeState(1)
  alive && add_part!(X, :Life, live=1)
  for _ in 1:n
    v = add_part!(X, :V)
    add_part!(X, :Life, live=v)
    add_edge!(X, v, 1)
  end
  homomorphism(alive ? LiveCell : DeadCell, X; initial=(V=[1],))
end;

view_life(codom(living_neighbors(3; alive=false)))

# ## Create model by defining update rules

# A cell dies due to underpopulation if it has < 2 living neighbors
underpop = 
  TickRule(:Underpop, to_life, id(DeadCell); ac=[NAC(living_neighbors(2))]);

# A cell dies due to overpopulation if it has > 3 living neighbors
overpop = 
  TickRule(:Overpop, to_life, id(DeadCell); ac=[PAC(living_neighbors(4))]);

# A cell is born iff it has three living neighbors
birth = TickRule(:Birth, id(DeadCell), to_life; 
                 ac=[PAC(living_neighbors(3; alive=false)),
                     NAC(living_neighbors(4; alive=false)),
                     NAC(to_life)]); # this rule does not apply if cell is alive

GoL = ABM([underpop, overpop, birth]);  # ABM is constituted by its transition rules
GoL_coords = AddCoords(GoL); # migrate ABM into schema with coordinates

# Create an initial state
G = make_grid([0 1 0;
               1 1 1;
               0 1 0])
view_life(G)

# Let's check that our rules work the appropriate way
#
# There are 12 dead cells and 13 live ones. 

match_coords(f::ACSetTransformation) =  codom(f)[f[:V](1), :coords]
match_coords(rule::ABMRule, X) = match_coords.(get_matches(AddCoords(rule), X))

match_coords.(homomorphisms(LiveCell|>AddCoords, G))

# Let's calculate how many matches we
# have for the rule which determines which cells die from underpopulation: 
match_coords(underpop, G)

# This is right, the corners have zero living neighbors and the bottom middle 
# cell only has one.
# 
# Now, let's calculate how many matches we
# have for the rule which determines which cells die from overpopulation: 

match_coords(overpop, G)

# Below are the coordinates of the 2 dead cells that will come to life:

match_coords(birth, G)

# Run the ABM
res = run!(GoL_coords, G);

# View resuls

imgs = view(res, view_life);

# We can see our starting point 

imgs[1]

# The next step 

imgs[2]

# The next step 

imgs[3]

# 

imgs[4]

# 

imgs[5]

# 

imgs[6]

# 

imgs[7]
# 

imgs[8]
# 

imgs[9]
# 

imgs[10]
# 

imgs[11]
# 

imgs[12]
# 

imgs[13]
# 

imgs[14]

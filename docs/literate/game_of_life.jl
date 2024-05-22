# # Conway's Game of Life
# ## Set-up
#
# The first step of running a Julia program is to load the external libraries 
# one will be using. We do this with a `using` statement.

using AlgebraicABMs, Catlab, AlgebraicRewriting

# We will turn on debug messages for AlgebraicABMs, which means we get to see 
# what the ABM is thinking at each time step when we eventually run it.
ENV["JULIA_DEBUG"] = "AlgebraicABMs"; 

#=
## Schema 

Defining an schema is stating what data is required to specify a state of the
simulation at some point in time. In AlgebraicJulia, this is done via declaring
a `Presentation`, i.e. a database schema. Objects (`Ob`, or tables) are types of
entities. Homs (`Hom`, or foreign keys) are functional relationships between the
aforementioned entities. AttrTypes are placeholders for Julia types, which are
assigned to objects via attributes (`Attr`).

The schema below *extends* the schema for symmetric graphs, which consists in
two tables (`E` and `V`, for edges and vertices) and homs (`src`, `tgt`) which
relate the edges to the vertices. We consider each vertex to be a cell in the
Game of Life. This is a generalization of the game, since it is not enforced
that every cell has eight neighbors. However the example we run at the end will
be a regular grid of vertices. 

Everything in the `@present` presentation below simply adds to the schema of
symmetric graphs, indicated by the subtype operator `<: SchSymmetricGraph`. We
need one more piece of information to specify a state of the world in the game
of life: which cells are alive? There are a few ways one could distinguish the
living cells from the dead ones. Here, we add a new table, `Life` which can be
though of as a set of "life tokens". A function, `live`, assigns each such token
to a cell, designating it to be alive. Thus the subset of living cells is the
image of the `live` function. 

To summarize, a state of the world in the Game of Life is a set of cells (i.e.
vertices), a set of edges (with `src`/`tgt` functions) which encode which cells 
are near each other, and a distinguished subset of cells given by the `live`
function which marks which cells are alive. This information is summarized in 
the following graphical depiction of the schema.
=#

@present SchLifeGraph <: SchSymmetricGraph begin 
  Life::Ob
  live::Hom(Life,V)
end

to_graphviz(SchLifeGraph); # visualize the schema


#=
If `SchLifeGraph` is the piece of data which describes what it means to be a 
world state, we need a Julia datatype whose values are world staets.

The `@acset_type` function takes the specification of a datatype from a schema, 
and it generates an efficient Julia datatype, which is essentially an in-memory
database with the schema given by `SchLifeGraph`. 

This line defines `LifeState` to be this new Julia type, and futhermore it
states that it satisfies the `AbstractSymmetricGraph` interface: Catlab already 
has a lot of generic code to work with SymmetricGraph-like things, so now we are 
free to use it with LifeState.
=#

@acset_type LifeState(SchLifeGraph) <: AbstractSymmetricGraph;

#=
Before we start defining the dynamics of our Game of Life, it'll be helpful to
visualize states of the world. One thing that `SchLifeGraph` did not have was  
(x,y) coordinate information. Because we were representing proximity of vertices
by the prescence or absence of an edge connecting them, these coordinates are
irrelevant to the model dynamics. 

However, the coordinates are useful for visualization, so the below code creates
a new schema which extends the previous one to add coordinates.
=#
@present SchLifeCoords <: SchLifeGraph begin
  Coords::AttrType
  coords::Attr(V, Coords)
end

to_graphviz(SchLifeCoords)

#=
Likewise for this new schema, we need to use `@acset_type` to turn the
description of the required data of a world state (the schema) into an actual
Julia datatype. One difference here is that there is now an attribute present,
which means we need to indicate what Julia type it should be. 

Our new Julia type is `LifeStateCoords`, and we pick `Tuple{Int,Int}` as our
representation of (x,y) coordinates.
=#
@acset_type LifeStateCoords(SchLifeCoords){Tuple{Int,Int}} <: AbstractSymmetricGraph;

#= 
The following code is helpful for visualizing a game state. It's not important 
to understand this code in order to understand the model. It uses coordinates 
if the input is a `LifeStateCoords`, but otherwise it places the vertices in 
an arbitrary location.
=#
function view_life(X::Union{LifeState, LifeStateCoords}, pth=tempname())
  pg = PropertyGraph{Any}(; prog="neato", graph=Dict(),
    node=Dict(:shape => "circle", :style => "filled", :margin => "0"),
    edge=Dict(:dir => "none", :minlen => "1"))
  add_vertices!(pg, nparts(X, :V))
  for v in vertices(X)
    is_alive = isempty(incident(X, v, :live))
    set_vprop!(pg, v, :fillcolor, is_alive ? "red" : "green")
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
  return G
end;

#=
This helper function creates an ordinary, square grid and initializes cells as
alive or dead via a boolean-valued input matrix. Like above, one does not need
to understand this in order to understand the Game of Life model - it's just a 
convenient way to generate some initial states to run the model on.
=#
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

#= 
Another helper function, which creates a game state on a square n × n grid.
If the keyword `random` is `true` (which it is, by default), then cells have a
50% chance of being alive.
=#
make_grid(n::Int, random=true) = 
  make_grid((random ? rand : zeros)(Bool, (n, n)));

view_life(make_grid(3)) # For example, visualize a random 3x3 grid

#=
We've now constructed two schemas and two datatypes. No doubt they are related
to each other, and AlgebraicJulia gives us the ability to take advantage of that
relationship in order to automatically values of `LifeStateCoords` into values
of `LifeState` (it does this is the obvious way: it throws away the coordinate
information). More interestingly, it gives us a way of converting a `LifeState`
into a `LifeStateCoords`: we give "generic" coordinates to each of the vertices.
=#

AddCoords = Migrate(SchLifeGraph, LifeState, SchLifeCoords, LifeStateCoords; 
                    delta=false);

# ## Building blocks for the Game of Life model 

#= 
We are now ready to build our model, which will have three rules: 
underpopulation, overpopulation, and birth. When we write these rules down, it 
will be very succinct due to the hard work of this section: laying out the 
building blocks. 

Our first building block is an instance of `LifeState`. Now, 
we don't want to think of this as a particular state of the game of life, but 
rather we want to think of it as a *pattern* for which we could look for 
*matches* inside a real game state.

This pattern is very simple: it consists of a single vertex. A *match* for 
this pattern in some world state `X` is just a choice of a vertex in `X`.
=#
const Cell = LifeState(1) # this syntax means: one vertex, nothing else

#=
Our next building block is also a pattern. It is a cell but, furthermore, it is 
a cell which has been marked as alive. This means it consists in one vertex,
one 'life token', and the function `live` which maps `Life₁` to `V₁`.
=# 
const LiveCell = @acset LifeState begin V=1; Life=1; live=1 end # a single living cell
const to_life = homomorphism(Cell, LiveCell)  # the unique map Dead → Live

# `to_life` is an important morphism to understand. 

PAC(m) = AppCond(m; monic=true) # Positive Application condition
NAC(m) = AppCond(m, false; monic=true) # Negative Application condition
TickRule(name, I_L, I_R; ac) = # Rule which fires on 1.0, 2.0, ...
  ABMRule(name, Rule(I_L, I_R; ac), DiscreteHazard(1); schema=SchLifeGraph);

"""Create a context of n living neighbors for either a dead or alive cell"""
function living_neighbors(n::Int; alive=true)::ACSetTransformation
  X = LifeState(1)
  alive && add_part!(X, :Life, live=1)
  for _ in 1:n
    v = add_part!(X, :V)
    add_part!(X, :Life, live=v)
    add_edge!(X, v, 1)
  end
  homomorphism(alive ? LiveCell : Cell, X; initial=(V=[1],))
end;

view_life(codom(living_neighbors(3; alive=false)))

# ## Create model by defining update rules

# A cell dies due to underpopulation if it has < 2 living neighbors
underpop = 
  TickRule(:Underpop, to_life, id(Cell); ac=[NAC(living_neighbors(2))]);

# A cell dies due to overpopulation if it has > 3 living neighbors
overpop = 
  TickRule(:Overpop, to_life, id(Cell); ac=[PAC(living_neighbors(4))]);

# A cell is born iff it has three living neighbors
birth = TickRule(:Birth, id(Cell), to_life; 
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

match_coords.(homomorphisms(AddCoords(LiveCell), G))

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

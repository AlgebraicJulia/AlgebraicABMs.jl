  # # Labor Market Search and Matching
# ## Set-up
#
# First, we load the necessary libraries from AlgebraicJulia and elsewhere.

using AlgebraicABMs, Catlab, AlgebraicRewriting, Random, Test, Plots, DataFrames, DataMigrations
using Catlab: to_graphviz # hide
import Distributions: Exponential, LogNormal
using Pipe: @pipe


ENV["JULIA_DEBUG"] = "AlgebraicABMs"; # hide
Random.seed!(123); # hide

# ## Background
# Inspired by the search-and-matching ABM implementation in [del Rio-Chanona et. al. (2022)](https://royalsocietypublishing.org/doi/suppl/10.1098/rsif.2020.0898), 
# we demonstrate a very basic ABM of a [labor market seach-and-matching model](https://en.wikipedia.org/wiki/Search_and_matching_theory_(economics)).  
# A state-of-the-world in our model will be represented by an ACSet on a schema that we define, and the model itself
# will consist of a set of rewrite rules for ACSets on that schema, each paired with a probability distribution over the
# waiting time before that rewrite "fires" on each of the matches it "recognizes" in the current state of the world.
# Running the model generates a trajectory of ACSets and times of state transitions.
#
# Since the ACSet rewrite rules are defined using ACSet instances and the homomorphisms (or ACSet Transformations) between
# them, it will be useful to revise AlgebraicRewriting.jl functions for creating them.  In practice, the rules below
# are defined using identitiy morphisms on coproducts of representables and homomorphism search for unique homomorphisms
# between them. 
# 
# ## Schema
# We define our Schema "from scratch" by specifying the types of objects in our model and the mappings 
# (or "homomorphisms") between them.  

@present SchLaborMarket(FreeSchema) begin
	Person::Ob
	Job::Ob
	Firm::Ob
	Vacancy::Ob

  employee::Hom(Job, Person)
  employer::Hom(Job, Firm)
  advertised_by::Hom(Vacancy, Firm)
end;

to_graphviz(SchLaborMarket) # hide


# We then create a Julia Type `LaborMarket` for instances of this schema.
@acset_type LaborMarket(SchLaborMarket);

# ## Constructing Instances
# Having defined the schema, we will build our model(s) by constructing particular 
# instances of this schema, and the transformations between them.  This can be done
# by figuring out how our desired instance would be constructed in memory, and adding
# the parts "by hand" using the imperative interface provided by `add_part` and related 
# functions.  As an alternative, we can specify some basic elements of our instances by 
# taking the freely constucted minimal example of each of our entities ("objects" - but 
# not in the sense of Object-Oriented Programming).  This creates a "generic" instance 
# of the chosen entity, which in particular doesn't force any two objects to be the same 
# when they don't have to be.
#  
# The "representable" Person and Firm are what we would expect - single instances of
# the relevant entity, and nothing else. 
P = representable(LaborMarket, :Person);
P |> elements |> to_graphviz # hide
#
F = representable(LaborMarket, :Firm);
F |> elements |> to_graphviz # hide


# The representable vacancy, however, can't be just a single vacancy and nothing else.
# To make a well-defined ACSet conforming to our schema, the representable Vacancy has 
# to have both a vacancy and a firm, with a function mapping the former to the latter.
V = representable(LaborMarket, :Vacancy);
V |> elements |> to_graphviz # hide

# The representable Job, in turn, has two functions pointing from it, so it has to
# include both a Person and a Firm.
J = representable(LaborMarket, :Job);
J |> elements |> to_graphviz # hide

# Joining these instances together by placing them "side by side" (i.e. not forcing any 
# entities from different ACSets to be equal to each other in the result) is known as 
# taking their coproduct, and has been implemented using the "oplus" symbol - $\oplus$.

generic_person_sidebyside_generic_firm = P⊕F
generic_person_sidebyside_generic_firm |> elements |> to_graphviz #hide
#  
one_of_each_generic_thing = P⊕F⊕V⊕J
one_of_each_generic_thing |> elements |> to_graphviz #hide

# The coproduct has a "unit" consisting of the "empty" instance of that schema.  
# We follow convention by denoting it with the letter O.

O = LaborMarket();
O |> elements |> to_graphviz # hide

# Instances that can be formed using oplus and the generic members of the objects are 
# known as "coproducts of representables".  One advantage of constructing our instances
# in this way is that we always know we're dealing with well-formed ACSets, which prevents
# cryptic errors further down the line.  However we may want to be able to express
# situations where the same entity plays more than one role in an instance.  We can still
# do this by constructing a free ACSet on a number of generic objects subject to equality
# constraints (known as a "colimit of representables").  The macro `@acset_colim` allows
# us to do this, if we give it a cached collection of all of the representables for our schema. 

yF = yoneda_cache(LaborMarket);


employer_also_hiring = @acset_colim yF begin
	j1::Job
	v1::Vacancy 
	employer(j1) == advertised_by(v1)
end; 

employer_also_hiring |> elements |> to_graphviz # hide

# ## Rules
# Now that we are able to construct instances of our schema - which we can think of as
# states of (part of) the world at given points in time - we can define the types of 
# change which can occur in our model.  These will take the form of ACSet rewriting rules,
# which are a generalization of graph rewriting rules (since a graph can be defined as a
# relatively simple ACSet, or indeed CSet).  We will use Double Pushout and Single Pushout
# rewriting, which both take an input pattern of the form  L ↢ I → R , where L is the input 
# pattern to be matched in the existing state of the world ("Before"), R is the output
# pattern that should exist going forward ("After") and I is the pattern of items in the
# input match which should carry over into the output match.
#
# The arrows in this pattern represent "ACSet transformations", which are structure-preserving
# morphisms between ACSets. The rewrite rule is specified using a pair of ACSet transformations
# (I $\rightarrowtail$ L and I → R) of ACSets sharing the same schema, where both transformations 
# have the same ACSet as their domain.  While the ACSet Transformations can be built using the 
# machinery available in [Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl), we will find in
# many cases that the transformations (also known as homomorphisms) between
# two relatively simple instances will be unique, so we only need to specify the domain
# and codomain ACSets and rely on homomorphism search to find the mapping we intend.  In 
# the case where the domain and codomain are the same, such as where the input pattern
# persists in its entirety, we can specify the "transformation" mapping everything to
# itself using id().  If the left homomorphism is id(), then we won't delete any part
# of the pattern, if the right homomorphism is id(), then we won't add anything that wasn't
# there before.

# One of the events that can occur in our model is that any firm which exists can post
# a vacancy.  The input pattern is the representable firm (L = F), the firm persists (I = F),
# and in the output pattern, it has become part of a connected vacancy-firm pair (R = V).

post_vacancy = Rule{:DPO}(
	id(F),
	homomorphism(F, V)
);

# Swapping the pattern has the effect of running the rule in reverse.

withdraw_vacancy = Rule{:DPO}(
	homomorphism(F, V),
	id(F)
);


# People can appear out of nothing ...
birth = Rule{:DPO}(
  id(O),
  homomorphism(O, P)
);

# ... and unto nothing they shall return.  This rule uses Single Pushout rewriting,
# which uses "cascading deletes" to eliminate entities which are mapped by functions
# to an entity which gets deleted.  We use SPO for this reule because we want to 
# eliminate any jobs which point to the now-defunct person.
death = Rule{:SPO}(
	homomorphism(O, P),
	id(O)
);

# Firms come and go in the same way, like Soho Italian restaurants in a Douglas Adams novel.

firm_entry = Rule{:DPO}(
  id(O),
  homomorphism(O, F)
);

# If we apply the firm entry rule to an empty world-state, the result is a world that has
# exactly one firm and nothing else.

@assert rewrite(firm_entry, O) == F

# And likewise, in a world state with a single firm, if that firm exits we are left with
# the empty world state.

firm_exit = Rule{:SPO}(
	homomorphism(O, F),
	id(O)
);

@assert rewrite(firm_exit, F) == O
 
# To make an ABM, we wrap a rule in a named container with a probability distribution over
# how long it takes to "fire".  A model is created from a list of these wrapped rules.  To
# demonstrate, we make a trivial model to show a population converging to a steady state based
# on constant birth and mortality hazards.  

birth_abm_rule = ABMRule(:Birth, birth, ContinuousHazard(1/50));
death_abm_rule = ABMRule(:Death, death, ContinuousHazard(1));

people_only_abm = ABM([birth_abm_rule, death_abm_rule]);

# 

function sequence_of_states(results::AlgebraicABMs.ABMs.Traj)             # hide
  cat([results.init], [codom(right(h)) for h in results.hist], dims=1)    # hide
end 																																			# hide

function event_times(results::AlgebraicABMs.ABMs.Traj)									  # hide
  cat([0.0], [e[1] for e in results.events]; dims=1)											# hide
end 																																			# hide

function unpack_results(r::AlgebraicABMs.ABMs.Traj) 											# hide
  DataFrame(																															# hide
  	time = event_times(r),																							  # hide	
  	state = sequence_of_states(r)																				  # hide	
  )                                                                       # hide
end; 																																			# hide

function obj_counts(abm_state::LaborMarket)																# hide
  [k => length(v)																													# hide
    for (k, v) in 																											  # hide		
    zip(keys(abm_state.parts), abm_state.parts)														# hide
  ] |> NamedTuple 																												# hide
end; 																																			# hide

function full_df(results::AlgebraicABMs.ABMs.Traj)												# hide
  state_time_df = unpack_results(results)																	# hide
  obj_counts_df = obj_counts.(state_time_df.state) |> DataFrame  					# hide
  hcat(state_time_df, obj_counts_df)																			# hide
end; 																																			# hide

function plot_full_df(df::DataFrame)																			# hide
  Plots.plot(																															# hide
  	df.time,																															# hide
  	[df.Person, df.Firm, df.Job, df.Vacancy];															# hide
  	 labels=["Person" "Firm" "Job" "Vacancy"]															# hide
  )																																				# hide
end; 																																			# hide

function plot_full_df(results::AlgebraicABMs.ABMs.Traj)										# hide
	plot_full_df(full_df(results))																					# hide
end; 																																			# hide

# If we run the same ABM on an initial state which has Firms in it, they just sit there, untouched by
# either of the ABM rules.


# To give the firms their own dynamics, we include two more ABM rules, using the firm entry and exit
# patterns defined above.  We then add two more rules to generate a steady state of vacancies.

people_and_firms_abm = @pipe people_only_abm |>
  copy(_) |>
  push!(_, ABMRule(:FirmEntry, firm_entry, ContinuousHazard(1/10))) |>
  push!(_, ABMRule(:FirmExit, firm_exit, ContinuousHazard(1))) |>
  push!(_, ABMRule(:PostVacancy, post_vacancy, ContinuousHazard(1))) |>
  push!(_, ABMRule(:WithdrawVacancy, withdraw_vacancy, ContinuousHazard(1)));


# Hiring, however, presents a new challenge.  We want to convert a
# person and vacancy-firm pair to a person and firm connected by a job, with the person
# and firm staying the same, and the vacancy disappearing.  But we don't want to keep 
# adding jobs to the same person indefinitely - in this case, we'll abstract from reality a little,
# and pretend that people can have at most one job.  We enforce this using an "application condition"
# attached to our rule.  This is formed by specifying a further homomorphism from L to another pattern,
# and specifying whether that further match is required or forbidden (false).  One could say it puts the 
# pattern in a larger context, and the boolean flag says whether this larger context is either mandatory 
# or forbidden.  In this case, we want to rule out situations where the Person and Firm-Vacancy pair are 
# part of a larger pattern consisting of a Job-Person-Firm triple and a Firm-Vacancy pair - i.e. situations
# in which the person we are matching already has a job.

hire = Rule{:DPO}(
	homomorphism(P⊕F, P⊕V),
	homomorphism(P⊕F, J);
	ac = [
	  AppCond(homomorphism(P⊕V, J⊕V), false) # Limit one job per person
	]
);


# Separations occur when we match a connected Person-Firm-Job triple and the person
# and firm continue separately but the job disappears. 
fire = Rule{:DPO}(
	homomorphism(P⊕F, J),
  id(P⊕F)
);

# We assume at first that hiring and firing occur at constant rates.
constant_job_dynamics_abm = @pipe people_and_firms_abm |>
  copy(_) |>
  push!(_, ABMRule(:Hire, hire, ContinuousHazard(1))) |>
  push!(_, ABMRule(:Fire, fire, ContinuousHazard(1)))

# In particular, we care about how many people don't have jobs.
function number_unemployed(LM::LaborMarket)
	nparts(LM,:Person) - nparts(LM,:Job)
end;

# Following the literature, we measure the interplay of supply and demand
# in the labour market using this "market tightness" ratio.

function market_tightness(state_of_world::LaborMarket)
  length(state_of_world.parts.Vacancy)/number_unemployed(state_of_world)
end;

# This may come in handy for defining distributions
function logistic_function(L, k)
  (x -> L / (1 + exp(-k*x)))
end;

# In order to make our model a little more interesting, we can make the likelihood of a worker
# filling a vacancy depend on supply and demand via a function of market tightness, as defined
# above.  We define a time-independent function of the match which returns a distribution over 
# firing times, to replace the (constant) function defined by ContinuousHazard.  We can then use 
# that as the firing distribution for our Hire rule.

dependent_match_function = (
	m -> begin
	  state_of_world = codom(m)
	  v_over_u = market_tightness(state_of_world)
	  q = logistic_function(2, 1)(v_over_u)  # decreasing function, elasticity between 0 and -1
	  Exponential(1/q)
  end
) |> ClosureState; 

full_abm = @pipe constant_job_dynamics_abm |>
  copy(_) |>
  filter(r -> r.name != :Hire, _) |> # Get rid of old hire rule
  push!(_, ABMRule(:Hire, hire, dependent_match_function)) # add new one

# ## Running the Model
# We can then construct an acset to reflect our starting state, and run a simulation for 
# a fixed amount of simulation time (as we do in this case) or until a fixed number of 
# events have happened.  NB - pending an update to the Single Pushout rewriting capability 
# ofAlgebraicABMs, we have temporarily deactivated the firm and person birth and death rules. 

partial_abm = filter(full_abm) do r 
	!(r.name in [:Birth, :Death, :FirmEntry, :FirmExit])
end;

initial_state = @acset LaborMarket begin
	Person = 20
	Firm = 6
end;

initial_state |> elements |> to_graphviz # hide

# 
result = run!(
  partial_abm,
  initial_state,
  maxtime = 100
);

plot_full_df(result) # hide

# 
function plot_beveridge_curve(results::AlgebraicABMs.ABMs.Traj) # hide
	states = sequence_of_states(results)                          # hide
  Plots.plot(                                                   # hide
    [number_unemployed(s)/nparts(s, :Person) for s in states],  # hide
    [market_tightness(s) for s in states],                      # hide
    xlabel = "U rate", ylabel = "V/U"														# hide
  )																														  # hide	
end 																														# hide

plot_beveridge_curve(result) # hide
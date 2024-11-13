# # Petri-net based ABMs
# 
# First we want to load our packages with `using

using AlgebraicABMs, Catlab, DataMigrations
using StateCharts
using Makie, CairoMakie


ENV["JULIA_DEBUG"] = "AlgebraicABMs";  # turn on @debug messages for this package

# SIR Agent-based model
# This model is the same  as the AnyLogic example model "SIR Agent-Based Networks.". 
# This model currently support three types of networks:
#   1. Random (if the parameter "p_random_connect" = 1.0)
#   2. Ring lattice (if the parameter "p_random_connect" = 0.0)
#   3. Small world (if the parameter "p_random_connect":  0.0 < p_random_connect < 1.0 )
#       it indicates that there are "p_random_connect" percentage connections are radom connections, and 
#       ( 1 - "p_random_connect" ) percentage connections are connections are Ring lattice connections.
# This example model in Anylogic can be found: https://cloud.anylogic.com/model/7088e817-1dab-42b1-89fe-b19bc0a823e1?mode=SETTINGS

# Step 1: define parameters
# parameters of state transition
total_population = 50 # Unit: persons
frac_initial_infected = 0.02
contact_rate = 5 # Unit: contacts per day
infectivity = 0.01
average_illness_duration = 15 # days
# parameters define the network
average_connections = 2
p_random_connect = 0.9

# Step 2: create the StateChart

# function UnitStateChartF takes in 3 input argummapents:
# states: tuple/array of all states
# transitions: tuple/array of all transitions. The syntax of each transition:
#              transition_name => (source_state => target_state, tranisiton_type => values)
# alternative transitions: tuple/array of all alternative transitions. The syntax of each altanative transition:
#              :source_transition_name=>((alternative_transition_name, probability_value) => target_state)
states = [:S, :I, :R] # define the states
transitions = [ :Infection => (:S => :I, :Pattern => 1), :Recovery => (:I => :R, :Rate => 1.0 / average_illness_duration)]
alternatives = [] # this model does not include alternative transitions
# 2.1 create the Statechart
SIRStatechart = UnitStateChartF(states, transitions, alternatives)
# 2.2 Visualization of the Statechart
stateColors = Dict(:S => "green",:I => "red",:R => "gray"); # argument define colors of each state in StateChart
StateCharts.Graph(SIRStatechart, stateColors = stateColors)

# Step 3: create the model schema
### generate the model schema by composing auto state chart schema and the network schema

# 3.1 define the persons object named :V, because we plan to compose the schema_statechart with graph schema by identifying ":V"
schema_statechart = StateChartABMSchema_MultipleObjects(SIRStatechart,:V) |> schemaACSet |> schemaPresent
# 3.2 compose the state chart schema with the network schema (symmetric reflective graph) as the ABM model schema (without equations)
schema_model_without_equations = compose(Open([:S],schemaACSet(schema_statechart),[:V]), Open([:V],schemaACSet(SchUndirectedReflectiveNetwork),[:E])) |> apex |> schemaPresent 
# 3.3 add the composition equations of the model schema, since those equations disapper after composition. This would be fixed in the future using GATlab
@present schema_model <: schema_model_without_equations begin
    compose(inv,inv) == id(E)
    compose(inv,src) == tgt
    compose(inv,tgt) == src
    compose(refl, src) == id(V)
    compose(refl, tgt) == id(V)
    compose(refl, inv) == refl 
end

@acset_type SIRModelStaticNet(schema_model, index=[:src,:tgt]){Symbol} <: AbstractSymmetricReflexiveGraph
# show the model schema
to_graphviz(schema_model; prog="dot")

# Step 4: Data Migration from the ACSet with schema -- schema_statechart to the ACSet with schema -- schema_model
# 4.1 define the data migration rule to automatically migrate ACSets from pure state chart schema to the composed model shema    
const migrate_rule = @migration schema_model schema_statechart begin
    V => V; 
    ID=>ID; VID => VID
    S => S; SV => SV
    I => I; IV => IV
    R => R; RV => RV

    E => @product begin
        p1::V
        p2::V
    end
    src => p1
    tgt => p2
    inv => begin
        p1 => p2
        p2 => p1
    end
end

# 4.2 define the user_defined rewrite rules: rewrite rules for contacts of Infectives
# Note that each transitions_rule has a pair: timer=>rule
transition_rules = [ ContinuousHazard( 1.0 / (contact_rate * infectivity)) => make_infectious_rule_MultipleObjects(SIRStatechart, [[:S],[:I]],[:I],[[:I],[:I]], :V; use_DataMigration = true, acset = SIRModelStaticNet, migration_rule = migrate_rule)]

# Step 5: define the network
network = smallworldNetWork(Int(total_population), average_connections, p_random_connect);

# Step 6: Initialization
# In the initialization, randomly assign "total_population * frac_initial_infected" persons as Infectives (state I), and the rest persons are all Susceptibles (state S)
init = radomlyAssignInitialInfectives(SIRModelStaticNet(), Int(total_population), Int(frac_initial_infected * total_population); obn=:V, network = network)

# Step 7: run the ABM model
ts = 50.0 # the total running time
abm = make_ABM(SIRStatechart, transition_rules, :V; is_schema_singObject=false, use_DataMigration=true, acset = SIRModelStaticNet, migration_rule = migrate_rule)
res = run!(abm, init; maxtime=ts);

# Step 8: Visualization of results
# 8.1 plot out the time series of each state
Makie.plot(res; Dict(o=>X->nparts(X,o) for o in states)...)
# 8.2 show the networks at time t
t0 = 0 # initial state
StateCharts.Graph(res,t0,stateColors)

t10 = 10 # initial state
StateCharts.Graph(res,t10,stateColors)

t20 = 20
StateCharts.Graph(res,t20,stateColors)

t30 = 30
StateCharts.Graph(res,t30,stateColors)

t40 = 40
StateCharts.Graph(res,t40,stateColors)

t50 = Int(ts) # final state
StateCharts.Graph(res,t50,stateColors)
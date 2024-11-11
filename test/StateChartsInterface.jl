module TestStateChartsInterface

using Test
using Makie, CairoMakie
using StateCharts
using AlgebraicABMs
using AlgebraicABMs.ABMs: RepresentableP, RegularP
using Catlab


totalPopulation=5.0
β = 0.1
p_becomingInfective = 5.0

# create the pertussis state chart
pertussis = UnitStateChartF((:S,:E,:I,:R₄,:V₄), # states
(:newExposure=>(:S=>:E,:Pattern=>1),:becomingInfective=>(:E=>:I,:TimeOut=>p_becomingInfective),:recovery=>(:I=>:R₄,:Rate=>1.0/21.0),:vaccinated=>(:S=>:V₄,:Rate=>0.01)), #non-start transitions
() # alternatives for non-start transitions
)

StateCharts.Graph(pertussis)
schema = StateChartABMSchema_MultipleObjects(pertussis) |> schemaACSet |> schemaPresent
to_graphviz(schema; prog="dot")


########### use ABM model schema of one object P #####################
########### the states of state charts are attributes ################

# define the new_infectious rule
transition_rules = [ ContinuousHazard( 1.0 / β) => make_infectious_rule_MultipleObjects(pertussis, [[:S],[:I]],[:I],[[:E],[:I]], RdoubleI = false)]

# Initial state for schema of multiple object
@acset_type PertussisNet(schema, index=[:src,:tgt]){Symbol}
init = radomlyAssignInitialInfectives(PertussisNet(), Int(totalPopulation), Int(1))

## test ABM rules
abm = make_ABM(pertussis, transition_rules; is_schema_singObject=false)

@test abm.rules[1].pattern_type isa RepresentableP
@test abm.rules[2].pattern_type == RegularP()
@test abm.rules[3].pattern_type isa RepresentableP
@test abm.rules[4].pattern_type isa RepresentableP
# transition "newExposure"
@test nparts.(Ref(codom(Catlab.left(abm.rules[1].rule))), [:S,:E,:I,:R₄,:V₄,:P]) == [1,0,1,0,0,2] # S->p1, I->p2
@test nparts.(Ref(dom(Catlab.left(abm.rules[1].rule))), [:S,:E,:I,:R₄,:V₄,:P]) == [0,0,1,0,0,2] # p1, I->p2
@test nparts.(Ref(codom(Catlab.right(abm.rules[1].rule))), [:S,:E,:I,:R₄,:V₄,:P]) == [0,1,1,0,0,2] # E->p1, I->p2
# transition "recovery"
@test nparts.(Ref(codom(Catlab.left(abm.rules[3].rule))), [:S,:E,:I,:R₄,:V₄,:P]) == [0,0,1,0,0,1] # I->p
@test nparts.(Ref(dom(Catlab.left(abm.rules[3].rule))), [:S,:E,:I,:R₄,:V₄,:P]) == [0,0,0,0,0,1] # p
@test nparts.(Ref(codom(Catlab.right(abm.rules[3].rule))), [:S,:E,:I,:R₄,:V₄,:P]) == [0,0,0,1,0,1] # R4->p
# transition "vaccinated"
@test nparts.(Ref(codom(Catlab.left(abm.rules[4].rule))), [:S,:E,:I,:R₄,:V₄,:P]) == [1,0,0,0,0,1] # S->p
@test nparts.(Ref(dom(Catlab.left(abm.rules[4].rule))), [:S,:E,:I,:R₄,:V₄,:P]) == [0,0,0,0,0,1] # p
@test nparts.(Ref(codom(Catlab.right(abm.rules[4].rule))), [:S,:E,:I,:R₄,:V₄,:P]) == [0,0,0,0,1,1] # V4->p

#res = run!(abm, init; maxtime=50.0)
#Makie.plot(res; Dict(o=>X->nparts(X,o) for o in [:S,:E,:I,:R₄,:V₄])...)

end
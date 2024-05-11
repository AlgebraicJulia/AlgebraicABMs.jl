module TestPetriInterface 

using Test
using AlgebraicABMs, Catlab, AlgebraicPetri

sir_pn= @acset LabelledPetriNet begin
  S=3; sname=[:S,:I,:R]
  T=7; tname=[:inf,:rec,:birth,:deathS,:deathI,:deathR,:wane]
  I=7; it=[1,1,2,4,5,6,7]; is=[1,2,2,1,2,3,3]
  O=5; ot=[1,1,2,3,7]; os=[2,2,3,1,1]
end
t1, t2, t3, t4, t5, t6, t7 = AlgebraicABMs.PetriInterface.make_rules(sir_pn)
@test nparts.(Ref(codom(left(t1))), [:S,:I,:R]) == [1,1,0] # S, I


end # module

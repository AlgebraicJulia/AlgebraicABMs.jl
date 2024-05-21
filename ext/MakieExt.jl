module MakieExt 

using AlgebraicABMs.ABMs: Traj
using Catlab: codom, right
import Makie: plot
using Makie: Figure, Legend, Axis, lines!

"""
Take a trajectory and draw lines based on functions ACSet -> Number
"""
function plot(t::Traj; kw...)
  fig = Figure();
  ax = Axis(fig[1,1])
  times = [0; first.(t.events)]
  states = [t.init; codom.(right.(t.hist))]
  Legend(fig[1, 2], map(collect(values(kw))) do f
    lines!(ax, times, f.(states))
  end, string.(collect(keys(kw))))
  return fig
end 

end # module
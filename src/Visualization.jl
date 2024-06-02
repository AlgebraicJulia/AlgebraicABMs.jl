module Visualization

import AlgebraicRewriting: view_traj 
using Catlab: codom, right

using ..ABMs
using ..ABMs: Traj

"""
Take a function which produces a graphviz representation. Makes 
a series of images.
"""
function Base.view(t::Traj, viewer; dirname="default")
  pth = mkpath(joinpath("traj",dirname)) # make the folder
  N = length(string(length(t)))
  fi(i::Int) = joinpath(pth, "$(lpad(i, N, "0")).svg")
  rm.(joinpath.(pth, readdir(pth))) # clear the folder
  states = [t.init; codom.(right.(t.hist))]
  data = [(0, "(Initial state)", nothing); t.events]
  map(enumerate(zip(states, data))) do (i, (state, (time, event, _)))
    i -= 1
    G = viewer(state, fi(i))
    G.graph_attrs[:label] = "Step $i @ t=$time: Event $event"
    G.graph_attrs[:labelloc] = "t"
    open(fi(i), "w") do io
      show(io, "image/svg+xml", G)
    end
    G
  end
end 

end # module

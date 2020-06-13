module MarkovModelAndMMs

include(joinpath(@__DIR__, "MarkovStuff.jl"))
using .MarkovStuff

include(joinpath(@__DIR__, "MMs.jl"))
using .MMs

include(joinpath(@__DIR__, "CostFunctions.jl"))
using .CostFunctions

export MarkovStuff, MMs, CostFunctions

end

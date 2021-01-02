module SimpleNavierStokes
export LidDrivenCavity

using LinearAlgebra, SparseArrays, IterativeSolvers

include("functions.jl")
include("lid_driven_cavity.jl")
include("plot_recipes.jl")

end # module

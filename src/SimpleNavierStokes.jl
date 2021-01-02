module SimpleNavierStokes
export LidDrivenCavity, ShowStreamlines

using LinearAlgebra, SparseArrays, IterativeSolvers, Plots

include("functions.jl")
include("lid_driven_cavity.jl")
include("plot_utils.jl")

end # module

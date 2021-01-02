module SimpleNavierStokes
export LidDrivenCavity

using LinearAlgebra, SparseArrays, IterativeSolvers

include("functions.jl")
include("lid_driven_cavity.jl")

struct Results
    ψ::Array
    ω::Array
    hist::Array
    x::Array
    y::Array
    tfinal
    steps
    Re
end

end # module

module SimpleNavierStokes
export LidDrivenCavity

using LinearAlgebra, SparseArrays

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



function LidDrivenCavity(;
                         tfinal = Inf,
                         Lx = 1, Ly = 1, CFL = 0.5, Re = 100,
                         Nx = 65, Ny = 65,
                         u_n = 1, u_s = 0, v_w = 0, v_e = 0,
                         printfreq = 10)
    t0 = time()
    println("------------------Ny = $(Ny), Nx = $(Nx) ---------------")
    Δy  = Ly/(Ny-1)
    Δx  = Lx/(Nx-1)
    x = 0:Δx:Lx
    y = 0:Δy:Ly
    Δt = CFL * Δx

    # Construct matrix for Poisson equation
    A_poisson = BuildPoissonMatrix(Ny,Nx,Δx,Δy) # for coNxgrad
    # Construct matrix for advection-diffusion equation
    ap,an,as,ae,aw = BuildAdvectionDiffusionCoefficients(Re,Δt,Δx,Δy)
    # allocate empty matrices for Gauss-Siedel solver
    Rp = zeros(Ny,Nx); res = zeros(Ny,Nx)

    # initialize ω and ψ
    ω = zeros(Ny,Nx)
    ψ = zeros(Ny,Nx)

    # keep track of changes with Δ's
    ω_old = zeros(Ny,Nx)
    ψ_old = zeros(Ny,Nx)
    ω_hist = []
    ψ_hist = []
    residual = 1

    ###########################
    k0,t = 0,0
    while t < tfinal && maximum(residual) > 1e-8
        t += Δt
        k0 += 1
        # Solve Poisson equation for ψ:
        LinearSolve!(A_poisson,ψ,-ω)

        VorticityBoundaryConditions!(ω,ψ,Δx,Δy,u_n,u_s,v_e,v_w)
        # Modify the explicit part of advection-diffusion equation
        BuildAdvectionDiffusionRHS!(Rp,ω,ψ,Δt,Δx,Δy,Ny,Nx,Re)

        # Solve advection-diffusion equation for ω:
        GaussSiedel!(ω,ap,an,as,ae,aw,Rp,res)

        residual = RecordHistory!.([ω,ψ],[ω_old,ψ_old],[ω_hist,ψ_hist])


        if (k0 % printfreq == 0)
            println("Step: $k0 \t Time: $(round(t,digits=3))\t",
                    "|Δω|: $(round((residual[1]),digits=8)) \t",
                    "|Δψ|: $(round((residual[2]),digits=8)) \t")
        end
    end
    tt = round(time() - t0, digits=3)
    println("This took $tt seconds.")
    println("--------------------------------------------------------")
    Results(ψ,ω,hcat(ω_hist,ψ_hist),x,y,t,k0,Re)
end
end # module

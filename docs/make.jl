using Documenter
using SimpleNavierStokes

makedocs(
    sitename = "SimpleNavierStokes.jl Documentation",
    pages = [
	     "Index" => "index.md",
	     "Gauss-Siedel Solver" => "gausssiedel.md",
	     ],
    format = Documenter.HTML(),
    modules = [SimpleNavierStokes]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

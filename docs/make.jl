using Documenter
using SimpleNavierStokes

makedocs(
    sitename = "SimpleNavierStokes",
    format = Documenter.HTML(),
    modules = [SimpleNavierStokes]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

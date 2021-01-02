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

ShowStreamlines(sol::Results) = contour(sol.x,sol.y,reverse(reverse(sol.ψ,dims=1),dims=2),
          aspectratio=1,framestyle=:box,
          xlims=(sol.x[1],sol.x[end]),
          ylims=(sol.y[1],sol.y[end]),
          legend=:none,grid=:none)

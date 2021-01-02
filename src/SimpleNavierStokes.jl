module SimpleNavierStokes
export foo, bar

foo(x::T, y::T) where T <: Real = x + y - 5
bar(z::Float64) = foo(sqrt(z), z)

greet() = print("Hello World!")

end # module

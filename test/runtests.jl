using SimpleNavierStokes, Test

@testset "Bare-bones test" begin
    @test LidDrivenCavity(tfinal=2) == 2
end


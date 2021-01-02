using SimpleNavierStokes, Test

@testset "Lid-Driven cavity runs" begin
    @test LidDrivenCavity(tfinal=1).tfinal == 1
    @test LidDrivenCavity(tfinal=2).tfinal == 2
end


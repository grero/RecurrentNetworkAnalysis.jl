using Test
using RecurrentNetworkAnalysis


@testset "Basic" begin
    
    # a single rotational mode
    a = [1.0-0.5im, 1.0+0.5im]

    pp = RecurrentNetworkAnalysis.get_rotational_modes(a)
    @test pp == [(2,1)]

    # skew symmetric matrix representing rotation
    A = [1.0 -0.01;0.01 1.0]

    pp = RecurrentNetworkAnalysis.get_rotational_modes(A)
    @test pp == [(2,1)]

    nn = RecurrentNetworkAnalysis.non_normality(A)
    @show nn == 0.0
end
include("../src/LinearProgramming.jl")
using .LinearProgramming
using Test

atol = 1e-5
scaling_methods = [nothing, "arithmetic", "geometric"]

function test_problem(A, c, b, symbols, true_X, true_z; approx::Bool=true)
    for scaling in scaling_methods
        lp = LinearProgram(A, c, b, symbols)
        X, z = revised_simplex(lp, scaling_method=scaling)
        println(X)
        println(z)
        if approx
            @test all(x -> x == 1, isapprox.(X, true_X, atol=atol))
            @test isapprox(z, true_z, atol=atol)
        else
            @test X == true_X
            @test z == true_z
        end
    end
end

@testset "LinearProgramming.jl" begin
    
    # Problem 1 - Simple problem
    A = [3 4 1; 4 1 2; 5 -3 2; -1 2 1]
    b = [16 12 10 4]
    c = [2 5 3]
    symbols = ["<=", "<=", "<=", "<="]
    true_X = [1.8333333333333333, 2.333333333333333, 1.1666666666666667]
    true_z = 18.83333333333333
    test_problem(A, c, b, symbols, true_X, true_z)

    # Problem 2 - Mixed constraints problem
    A = [3 4 1; 4 1 2; 5 -3 2; -1 2 1]
    b = [12 20 24 18]
    c = [1 5 3]
    symbols = [">=", "<=", "<=", "<="]
    true_X = [2.4444444444444446, 10.222222222222221, 0.0]
    true_z = 53.55555555555555
    test_problem(A, c, b, symbols, true_X, true_z)

    # Problem 3 - Mixed constraints problem 2
    A = [1 1 0; 1 0 1; 0 1 1]
    b = [20 5 10]
    c = [1 -1 3]
    symbols = ["<=", "=", ">="]
    true_X = [0.0, 5.0, 5.0]
    true_z = 10.0
    test_problem(A, c, b, symbols, true_X, true_z)

    # Problem 4 - Infeasible problem
    A = [1 0; 0 1; 1 1]
    b = [6 6 11]
    c = [1 1]
    symbols = [">=", ">=", "<="]
    true_X = []
    true_z = "Infeasible"
    test_problem(A, c, b, symbols, true_X, true_z, approx=false)

    # Problem 5 - Unbounded problem
    A = [6 5; 10 20; 1 0]
    b = [60 150 8]
    c = [500 450]
    symbols = [">=", ">=", ">="]
    true_X = []
    true_z = "Unbounded"
    test_problem(A, c, b, symbols, true_X, true_z, approx=false)

    # Problem 6 - Unbounded problem 2
    A = [1 -1; 1 1]
    b = [1 2]
    c = [1 1]
    symbols = [">=", ">="]
    true_X = []
    true_z = "Unbounded"
    test_problem(A, c, b, symbols, true_X, true_z, approx=false)

end
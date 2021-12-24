# LinearProgramming.jl

LinearProgramming.jl is a package that allows the user to solve linear programming problems using Julia.

Linear programming problems require the optimization of a linear **objective function**, subject to a series of linear equality and inequality **constraints**.
The typical linear programming problem has the following matrix form:

$$max \ \textbf{c}^T \cdot \textbf{x}$$
$$s.t. \ A\textbf{x} \leq \textbf{b} $$

This package solves such problems using the **revised simplex algorithm**, a much more efficient version of the simplex algorithm. The main challenges of this algorithm are achieving good computational performance and numerical stability. 

Numerical instability often arises because the coefficients in the objective function or in the constraints are badly scaled. In order to mitigate this issue, this package allows the user to apply a scaling step if desired. In particular, two scaling methods are currently supported:
- Arithmetic mean scaling
- Geometric mean scaling

Ploskas & Samaras (2013) show evidence that these two scaling methods are effective both to improve numerical stability and to reduce the number of iterations required by the algorithm to converge.

As for computational efficiency, this package achieves a performance comparable to that of well-known Julia linear programming packages, such as the one provided by JuMP. 

## Example

In order to solve a custom problem, the following steps are necessary.

Firstly, the linear program must be created. To do so, the LinearProgram type is used, which takes the following parameters:
- *A*: matrix of coefficients of the constraints
- *c*: vector of coefficients of the objective function
- *b*: right-hand side of the constraints
- *symbols*: a vector of equality/inequality symbols of the constraints

Consider this simple linear program as an example:

$$max_x \ z = 2x_1 + 5x_2 + 3x_3$$
$$s.t. \ \ 3x_1 + 4x_2 + 1x_3 \geq 12$$
$$\ \ \ \ \ \ \ 4x_1 + 1x_2 + 2x_3 \leq 20$$
$$\ \ \ \ \ \ \ 5x_1 - 3x_2 + 2x_3 \leq 24$$
$$\ \ \ \ \ \ \ -1x_1 + 2x_2 - 1x_3 \leq 18$$

This is how you create an instance of this problem with this package:

    # Define problem in matrix form
    A = [3 4 1; 4 1 2; 5 -3 2; -1 2 1]
    c = [2 5 3]
    b = [12 20 24 18]
    symbols = [">=", "<=", "<=", "<="]

    # Create the linear program
    lp = LinearProgram(A, c, b, symbols)

Once the linear program is instantiated, it can be easily solved using the *revised_simplex* function. This function turns the problem into canonical form, it then applies scaling if required and finally solves the problem. No scaling is applied by default. The function returns a vector of optimal values of X and the value of the objective function at the optimum. The algorithm can either find an optimal solution or mark the problem as unbounded/infeasible.

    # Find solution
    X, z = revised_simplex(lp, scaling_method="geometric")

The package also allows the user to find the **dual** linear program given its **primal**. The dual has the same properties and can be solved as the primal linear program.

    # Find dual linear program
    dual = dual_program(lp)

Currently it is not possible to compute reduced costs and shadow prices using this package. This is a limitation and is left as a future improvement. 

## References

- Bradley, Hux & Magnanti (1977), *Applied Mathematical Programming*
- Ploskas & Samaras (2013), *The impact of scaling on simplex type algorithms*
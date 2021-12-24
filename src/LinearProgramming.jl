module LinearProgramming



using LinearAlgebra



export 
# Types
    LinearProgram,
# Functions
    revised_simplex,
    dual_program,
    to_canonical_form!



"""
Data type for the linear program
"""
mutable struct LinearProgram
    A::Matrix{Float64}                  # coefficients in the constraints
    c::Matrix{Float64}                  # coefficients in the objective function
    b::Matrix{Float64}                  # RHS
    B::Matrix{Float64}                  # initial basis matrix
    cB::Matrix{Float64}                 # coefficients of basic variables in the original objective function
    symbols::Vector{String}             # (in)equality signs
    maximize::Bool                      # whether problem has to be maximized or minimized
    n_vars::Int                         # number of variables in the original problem
    n_cons::Int                         # number of constraints in the original problem
    canonical_form::Bool                # problem is in canonical form
    scaling_factors::Vector{Float64}    # scaling factors of the objective function coefficients
    function LinearProgram(A::Matrix, c::Union{Vector, Matrix}, b::Union{Vector, Matrix}, symbols::Vector{String}, maximize::Bool=true)
        # Checks on c
        typeof(c) <: Vector && (c = reshape(c, 1, length(c)))   # if c is a vector turn it into a matrix
        !maximize && (c *= -1)                                  # if function has to be minimized, multiply everything by -1
        n_vars = size(c)[2]
        # Checks on A
        size(A)[2] != n_vars && throw(ArgumentError("Size of constraint matrix incompatible with coefficient vector"))
        # Checks on b
        (typeof(b) <: Vector || size(b) == (1, length(b))) && (b = reshape(b, length(b), 1))
        n_cons = size(b)[1]
        n_cons != size(A)[1] && throw(ArgumentError("Size of vector of value constraints b does not match size of constraint matrix A"))
        # Checks on symbols
        length(symbols) != n_cons && throw(ArgumentError("Length of symbols vector different than number of constraints"))
        any(x -> x ∉ ("≤", "≥", "=", "<=", ">="), symbols) && throw(ArgumentError("Unsupported inequality symbol"))
        # Initialize basis matrix B
        B = Matrix{Int}(I, n_cons, n_cons)
        # Initialize vector of coefficients of the basic variables cB
        cB = zeros(1, n_cons)
        return new(A, c, b, B, cB, symbols, maximize, n_vars, n_cons, false, [])
    end
end



"""
Turn linear program into canonical form
"""
function to_canonical_form!(lp::LinearProgram)
    lp.canonical_form && throw(ArgumentError("Program is already in canonical form"))
    for i in 1:lp.n_cons
        # if RHS is negative multiply everything by -1 and flip symbol
        if lp.b[i] < 0
            lp.A[i,:] *= -1
            lp.b[i] *= -1
            if lp.symbols[i] ∈ ("≤", "<=")
                lp.symbols[i] = "≥"
            elseif lp.symbols[i] ∈ ("≥", ">=")
                lp.symbols[i] = "≤"
            end
        end
        # Handle case of ≥ symbol
        if lp.symbols[i] ∈ ("≥", ">=")
            lp.c = hcat(lp.c, 0)                                        # add a zero in the coefficient for each surplus variable added
            lp.cB[i] = -BigInt(2)^64                                    # use big-M method 
            lp.A = hcat(lp.A, [j == i ? -1 : 0 for j in 1:lp.n_cons])   # add surplus variable
        elseif lp.symbols[i] == "="                               
            lp.cB[i] = -BigInt(2)^64                            
        end
    end
    lp.canonical_form = true
end



"""
Find dual linear program given primal
"""
function dual_program(lp::LinearProgram)
    A = Matrix(lp.A')
    c = reshape(lp.b, 1, size(lp.b)[1])
    b = lp.c
    lp.maximize ? (maximize = false) : (maximize = true)
    symbols = String[]
    for i in 1:size(b)[2]
        if i <= size(c)[2] && lp.symbols[i] ∈ ("≥", ">=")
            push!(symbols, "≤")
        else
            push!(symbols, "≥")
        end
    end
    return LinearProgram(A, c, b, symbols, maximize)
end



"""
Compute a scaling factor for each row and column in the constraint coefficient matrix A. 
Then scale each row in A and each element in the RHS b by the row scaling factor and each column 
in A and each element in the objective function c by the column scaling factor.
"""
function scaling!(lp::LinearProgram, method::String="arithmetic")
    method ∉ ("arithmetic", "geometric") && throw(ArgumentError("Unsupported scaling method"))
    if method == "arithmetic"
        arithmetic_mean_scaling!(lp)
    else
        geometric_mean_scaling!(lp)
    end
end



"""
Arithmetic mean scaling as presented in Ploskas & Samaras (2013)
"""
function arithmetic_mean_scaling!(lp::LinearProgram)
    m, n = size(lp.A)
    A = copy(lp.A)
    # Row scaling
    for i in 1:m
        row = A[i,:]
        nz = row[row .!= 0]                       # non-zero elements of row i
        count_row, sum_row = length(nz), sum(abs.(nz))
        if count_row != 0 && sum_row != 0
            r = count_row / sum_row               # row scaling factor
            lp.A[i,:] = lp.A[i,:] * r             # scale constraint function coefficient
            lp.b[i] = lp.b[i] * r                 # scale RHS
        end
    end
    # Column scaling
    for j in 1:n
        col = A[:,j]
        nz = col[col .!= 0]
        count_col, sum_col = length(nz), sum(abs.(nz))
        if count_col != 0 && sum_col != 0
            s = count_col / sum_col               # column scaling factor
            lp.A[:,j] = lp.A[:,j] * s             # scale constraint function coefficient
            lp.c[j] = lp.c[j] * s                 # scale objective function coefficient
            push!(lp.scaling_factors, s)
        end
    end
end



"""
Geometric mean scaling as presented in Ploskas & Samaras (2013)
"""
function geometric_mean_scaling!(lp::LinearProgram)
    m, n = size(lp.A)
    A = copy(lp.A)
    # Row scaling
    for i in 1:m
        row = A[i,:]
        values = abs.(row[row .!= 0])             # abs of non-zero elements of row i
        min_row, max_row = minimum(values), maximum(values)
        if min_row != 0 && max_row != 0
            r = 1 / sqrt(min_row * max_row)       # row scaling factor
            lp.A[i,:] = lp.A[i,:] * r             # scale constraint function coefficient
            lp.b[i] = lp.b[i] * r                 # scale RHS
        end
    end
    # Column scaling
    for j in 1:n
        col = A[:,j]
        values = abs.(col[col .!= 0])
        min_col, max_col = minimum(values), maximum(values)
        if min_col != 0 && max_col != 0
            s = 1 / sqrt(min_col * max_col)       # column scaling factor
            lp.A[:,j] = lp.A[:,j] * s             # scale constraint function coefficient
            lp.c[j] = lp.c[j] * s                 # scale objective function coefficient
            push!(lp.scaling_factors, s)
        end
    end
end



"""
Solve linear programming problem using the Revised Simplex Method 
with Big-M Method to handle artificial variables (i.e. cases with ≥ and = signs)
"""
function revised_simplex(lp::LinearProgram; scaling_method=nothing, ϵ=1e-8)
    # Apply scaling if desired
    if !isnothing(scaling_method)
        scaling!(lp, scaling_method)
    end
    # Turn problem into canonical form if not already so
    if !lp.canonical_form
        to_canonical_form!(lp)
    end
    # Read model parameters
    A, c, b, B, cB, n_vars = lp.A, lp.c, lp.b, lp.B, lp.cB, lp.n_vars
    # Mapping from non-basic to basic variables used to retrieve optimal solution at the end
    d = Dict()
    previous_d = copy(d)         # mapping at the previous iteration, needed to stop the loop if problem is infeasible
    # Compute coefficients for the non-basic variables
    B_inv = inv(B)
    b_hat = B_inv * b
    c_hat = (cB * B_inv * A) - c
    # Choose variable to introduce into the basis
    s = getindex(argmin(c_hat), 2)
    # Check if we are in the simple case without ≥ signs
    if cB == zeros(1, length(b))
        abs(c_hat[s]) < ϵ && (c_hat[s] = 0)             # allow for some precision error
    end
    i = 0
    # Solve linear program
    while c_hat[s] < 0                                  # check optimality criterion
        A_hat_s = B_inv * A[:,s]
        # Choose variable to drop from the basis using the minimum ratio rule
        ratios = b_hat ./ A_hat_s
        ratios[ratios .≤ 0] .= Inf16                    # replace negative number with infinity
        # Check if problem is unbounded
        all(x -> x ≤ 0, A_hat_s) && (return ([], "Unbounded"))
        # Choose variable to drop from the basis
        r = getindex(argmin(ratios), 1)
        for pair in d                                   # substitute variables in basis
            r == pair[2] && delete!(d, pair[1])
        end
        # Update mapping between non-basic and basic variables
        d[s] = r                                       
        # Stop if problem is infeasible 
        d == previous_d && (return ([], "Infeasible"))
        previous_d = copy(d)
        # Update basis by replacing basic variable in row r with new variable s
        B[:,r] = A[:,s]
        # Update the coefficients of the basic variables in the objective function
        cB[r] = c[s]
        # Update coefficients for the non-basic variables
        B_inv = inv(B)
        b_hat = B_inv * b
        c_hat = (cB * B_inv * A) - c
        # Choose variable to introduce into the basis
        s = getindex(argmin(c_hat), 2)
        # Allow for some precision error
        abs(c_hat[s]) < ϵ && (c_hat[s] = 0)
    end
    # Compute optimal values of the decision variables
    X = zeros(n_vars)
    for var in d
        !(var[1] > n_vars) && (X[var[1]] = b_hat[var[2]])
    end
    # Unscale X to retrieve solution to the original problm if scaling was applied
    if !isnothing(scaling_method)
        X = X .* lp.scaling_factors
    end
    # Compute optimal value of the objective function
    z = cB * b_hat
    return X, z[1]
end



end
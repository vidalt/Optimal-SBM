### --------------------------------------------------------- ###
### MINLP

function MINLP(dataset::Dataset; 
               time_limit::Float64=600.0,
               symm_break_const::Bool=true,
               symm_w_const::Bool=true,
               w_lower::Float64=1e-12,
            )::Tuple{SBM, Matrix{Int}, OptResults}

    A = dataset.A
    K = dataset.n_communities
    n = dataset.n
    m = dataset.m
    k = dataset.k

    # Upper bound on w
    w_upper = getWUpperBound(dataset)
    @show(w_upper);

    # Set time limit
    open("couenne.opt", "w") do f
        write(f, "time_limit $time_limit\n");
        write(f, "allowable_gap 1e-6\n");
        write(f, "allowable_fraction_gap 1e-6\n");
        write(f, "acceptable_tol 1e-6\n");
        write(f, "tol 1e-6\n");
        write(f, "feas_tolerance 1e-7\n");
        write(f, "integer_tolerance 1e-7\n");
        write(f, "constr_viol_tol 1e-7\n");
        write(f, "compl_inf_tol 1e-7\n"); 
        write(f, "acceptable_constr_viol_tol 1e-7\n");
        write(f, "acceptable_compl_inf_tol 1e-7\n");
    end

    # --- Optimization model -------------------------------------------------
    model = Model(() -> AmplNLWriter.Optimizer("couenne"))

    # Variables
    @variable(model, w_lower <= w[1:K,1:K] <= w_upper);
    @variable(model, z[1:n,1:K], Bin);

    # Objective
    @NLobjective(model, Min, 0.5 * sum(
            (( -A[i,j]*log(w[r,s]) + ((k[i]*k[j])/(2*m))*w[r,s] ) * z[i,r] * z[j,s] )
            for r=1:K, s=1:K, i=1:n, j=1:n if A[i,j] != 0) + 
            0.5 * sum(
            (( ((k[i]*k[j])/(2*m))*w[r,s] ) * z[i,r] * z[j,s] )
            for r=1:K, s=1:K, i=1:n, j=1:n if A[i,j] == 0));

    # Constraints
    @constraint(model, con[i = 1:n], sum(z[i,r] for r=1:K) == 1); # assignment constraint
    
    # Symmetry breaking in numbered clustering - Frank Plastria
    if symm_break_const == true
        # object 1 must be in cluster 1
        @constraint(model, z[1,1] <= 1);
        @constraint(model, z[1,1] >= 1);
        @constraint(model, symBreak[r=2:(K-1), j=r:n], sum(z[i,l] for l=1:(r-1), i=2:(j-1)) - sum(z[j,l] for l=1:r) <= j - 3 )
    end

    # Enforce symmetry on w
    if symm_w_const == true
        @constraint(model, symL[r=1:K, s=1:(r-1)], w[r,s] - w[s,r] <= 0.0)
        @constraint(model, symG[r=1:K, s=1:(r-1)], w[r,s] - w[s,r] >= 0.0)
    end
    # ------------------------------------------------------------------------

    # --- Solve optimization model -------------------------------------------
    start = time();
    println("Start solving")
    optimize!(model)
    status = string(termination_status(model))
    if status == "LOCALLY_SOLVED"
        status = "OPTIMAL"
    elseif status == "OTHER_LIMIT"
        status = "TIME_LIMIT"
    end
    solvetime = time() - start
    obj_lb = Float64(objective_bound(model));
    obj_ub = Float64(objective_value(model));
    obj_lb = -Inf
    if (status == "OPTIMAL")
        obj_lb = obj_ub
    end
    
    # ------------------------------------------------------------------------

    # --- Recover variable values --------------------------------------------
    w_opt = value.(w)
    z_opt = value.(z)
    z_opt = round.(z_opt[:,:])
    z_opt[map(v -> isnan(v) , z_opt)] .= 0 # replace NaN values with 0
    z_opt = Int.(z_opt)
    # ------------------------------------------------------------------------

    # SBM
    sbm = SBM(w_opt, "poisson") # assumes a poisson distribution

    # Build OptResults
    nodecount = nothing
    lazycount = nothing
    iterations = nothing
    println()
    println("Obj. LB = $obj_lb")
    println("Obj. UB = $obj_ub")
    println("Status: $status")
    println("Solve time: $solvetime")

    # Print results
    @show(w_opt)
    @show(z_opt)
    
    opt_results = OptResults(obj_lb, obj_ub, status, solvetime, iterations, nodecount, lazycount, nothing)
    displayResults(opt_results)
    return sbm, z_opt, opt_results
end


### --------------------------------------------------------- ###
### MILP

function MILP(dataset::Dataset;                          # dataset
              time_limit::Float64=600.0,                 # time limit (in seconds)
              symm_break_const::Bool=true,               # whether to add symmetry-breaking constraints in the model
              symm_w_const::Bool=true,                   # whether to enforce symmetry on w
              balanced_communities::Bool=false,          # whether to enforce balanced communities
              TOL::Float64=1e-6,                         # tolerance
              w_lower::Float64=1e-12,                    # lower bound on w
              n_breakpoints::Int64=50)                   # number of breakpoints
    """
    Solves the Mixed Integer Linear Programming formulation.

    Parameters
    ----------
    dataset : Dataset
        Dataset representing an observed graph.

    Returns
    -------
    sbm : SBM
        Stochastic Block Model
    z : Array{Int64,2}
        Matrix of assignments of nodes to groups
    opt_results : OptResults
        Results of the optimization process
    """

    # --- Retrieve instance data -------------------------------------------- 
    A = dataset.A
    K = dataset.n_communities
    n = dataset.n
    m = dataset.m
    k = dataset.k

    println("K = $K groups\nn = $n nodes\nm = $m edges");
    println("Breakpoints: $n_breakpoints")


    # Calculate constant term
    Cnst = getObjectiveConstant(dataset)
    
    # Upper bound on w
    rho = 0
    for i=1:n, j=1:n
        if A[i,j] != 0
            rho = maximum([rho, (A[i,j] / (k[i]*k[j]))])
        end
    end
    w_upper = 2 * m * rho
    @show(w_upper)
    # ------------------------------------------------------------------------ 

    # --- Initialize all break-points ---------------------------------------- 
    # Initialize wp
    wp = collect(range(w_lower, stop=w_upper, length=n_breakpoints))
    # ------------------------------------------------------------------------ 


    # --- Linear pieces coefficients------------------------------------------ 
    # Calculating slope and intercept of linear approximations
    # Linear pieces coefficients
    # Calculating slope and intercept of linear approximations
    a = Array{Array{Float64},2}(undef, n, n) # Slope of each line
    b = Array{Array{Float64},2}(undef, n, n) # Intercept of each line
    for i=1:n, j=1:n
        if A[i,j] != 0
            # slope: a = f'(xp) ( derivative evaluated at the point wp[p] )
            a[i,j] = -A[i,j]./wp .+ (k[i]*k[j])/(2*m)
            # intercept: b = f(xp) - f'(xp)*xp ( derive it from yp = m * xp + b -> b = yp - m * xp )
            b[i,j] = A[i,j] .* (1 .- log.(wp))
        else
            a[i,j] = [(k[i]*k[j])/(2*m)]
            b[i,j] = [0]
        end
    end
    # ------------------------------------------------------------------------ 


    # --- Defining bounds on variables --------------------------------------- 
    M_lower = zeros(n,n,K,K);
    M_upper = zeros(n,n,K,K);
    for i=1:n, j=1:n, r=1:K, s=1:K
        if A[i,j] == 0
            M_lower[i,j,r,s] = 0 # Lower bound
            M_upper[i,j,r,s] = ((k[i]*k[j])/(2*m))*w_upper  # Upper bound
        else
            # Lower bound
            M_lower[i,j,r,s] = (A[i,j] * (1 - log(A[i,j]) + log((k[i]*k[j])/(2*m))))

            # Upper bound
            M_upper[i,j,r,s] = maximum([-A[i,j]*log(w_lower) + ((k[i]*k[j])/(2*m))*w_lower,
                                        -A[i,j]*log(w_upper) + ((k[i]*k[j])/(2*m))*w_upper])
        end
    end
    M_LOWER = minimum(M_lower)
    M_UPPER = maximum(M_upper)
    # ------------------------------------------------------------------------ 


    # --- Optimization model ------------------------------------------------- 
    solver = CPLEX.Optimizer
    model = Model(solver)
    set_optimizer_attribute(model, "CPXPARAM_MIP_Strategy_CallbackReducedLP", 0) # Access node information of the original model in the callback: CPX_OFF
    set_optimizer_attribute(model, "CPXPARAM_ScreenOutput", 1)              # Switching ON the display: CPX_ON
    set_optimizer_attribute(model, "CPXPARAM_Read_DataCheck", 1)            # Print warnings: CPX_DATACHECK_WARN
    set_optimizer_attribute(model, "CPX_PARAM_THREADS", 1)                  # Number of threads
    set_optimizer_attribute(model, "CPX_PARAM_TILIM", time_limit) # Time limit for the solver
    set_optimizer_attribute(model, "CPX_PARAM_EPGAP", TOL)                  # Relative MIP gap tolerance
    set_optimizer_attribute(model, "CPX_PARAM_PRELINEAR", 0)                # CPX_OFF
    set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)                   # CPX_PREREDUCE_NOPRIMALORDUAL
    set_optimizer_attribute(model, "CPX_PARAM_EPINT", 1e-9)                 # Integrality tolerance
    set_optimizer_attribute(model, "CPX_PARAM_EPRHS", 1e-9)                 # Feasibility tolerance
    set_optimizer_attribute(model, "CPX_PARAM_EPOPT", 1e-9)                 # Optimality tolerance
    set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Reduce", 1)      # Primal and dual preprocessing reduction type: Only primal reductions
    
    # Variables
    @variable(model, w_lower <= w[1:K,1:K] <= w_upper);
    @variable(model, z[1:n,1:K], Bin)
    @variable(model, 0 <= y[1:n,1:n,1:K,1:K] <= 1);
    @variable(model, M_LOWER <= x[1:n,1:n,1:K,1:K] <= M_UPPER);

    # Objective
    @objective(model, Min, - Cnst + 0.5 * sum( x[i,j,r,s] for r=1:K, s=1:K, i=1:n, j=1:n) );
    
    # Constraints
    # PWL Constraints
    @constraint(model, bigM1[i=1:n, j=1:n, r=1:K, s=1:K, p=1:length(a[i,j])], x[i,j,r,s] >= a[i,j][p]*w[r,s] + b[i,j][p] - M_upper[i,j,r,s]*(1 - y[i,j,r,s] )); # Big M-1 constraint
    
    # Big M constraints
    @constraint(model, bigM2[i=1:n, j=1:n, r=1:K, s=1:K], x[i,j,r,s] <= M_upper[i,j,r,s] * y[i,j,r,s] ); # Big M-2 constraint
    @constraint(model, bigM3[i=1:n, j=1:n, r=1:K, s=1:K], x[i,j,r,s] >= M_lower[i,j,r,s] * y[i,j,r,s] ); # Big M-3 constraint
    
    # Indicator variables constraints
    @constraint(model, bigM4[i=1:n, j=1:n, r=1:K, s=1:K], z[i,r] - y[i,j,r,s] >= 0 ); # Big M-4 constraint
    @constraint(model, bigM5[i=1:n, j=1:n, r=1:K, s=1:K], z[j,s] - y[i,j,r,s] >= 0 ); # Big M-5 constraint
    @constraint(model, bigM6[i=1:n, j=1:n, r=1:K, s=1:K], 1 - z[i,r] - z[j,s] + y[i,j,r,s] >= 0 ); # Big M-6 constraint

    # Assignment constraints
    @constraint(model, assign[i=1:n], sum(z[i,r] for r=1:K) == 1); # assignment constraint

    # Balanced communities
    if balanced_communities == true
        @constraint(model, bal[r=1:K], sum(z[i,r] for i=1:n) == n/K)
    end

    # Symmetry breaking constraints
    if symm_break_const == true
        @constraint(model, conObj1, z[1,1] == 1); # object 1 must be in cluster 1
        @constraint(model, symBreak[r=2:(K-1), j=r:n], sum(z[i,l] for l=1:(r-1), i=2:(j-1)) - sum(z[j,l] for l=1:r) <= j - 3 )
    end

    # Enforce symmetry on w
    if symm_w_const == true
        @constraint(model, symL[r=1:K, s=1:(r-1)], w[r,s] - w[s,r] <= 0.0)
        @constraint(model, symG[r=1:K, s=1:(r-1)], w[r,s] - w[s,r] >= 0.0)
    end

    # Get number of constraints of type GenericAffExpr
    num_cons = sum([num_constraints(model, func, set) 
                    for (func, set) in list_of_constraint_types(model) 
                    if func == GenericAffExpr{Float64, VariableRef}])
    println("Number of constraints before start: $num_cons")

    # Get the number of constraints of type VariableRef
    num_cons_vars = sum([num_constraints(model, func, set) 
                         for (func, set) in list_of_constraint_types(model) 
                         if func == VariableRef])
    println("Number of constraints on variables before start: $num_cons_vars")
    # ------------------------------------------------------------------------

    # --- Lazy constraint callback ------------------------------------------- 
    # Lazy constraints counter
    lazycount = 0

    # Defining lazy constraint callback
    function lazy_callback_function(cb_data)
        x_val = zeros(size(x))
        y_val = zeros(size(y))
        w_val = zeros(size(w))
        z_val = zeros(size(z))
        for i=1:n, j=1:n, r=1:K, s=1:K 
            x_val[i,j,r,s] = callback_value(cb_data, x[i,j,r,s])
        end
        for i=1:n, j=1:n, r=1:K, s=1:K 
            y_val[i,j,r,s] = callback_value(cb_data, y[i,j,r,s])
        end
        for r=1:K, s=1:K
            w_val[r,s] = callback_value(cb_data, w[r,s])
            if (w_val[r,s] <= 0)
                w_val[r,s] = 1e-50; # TODO: use a tolerance value here
            end
        end
        for i=1:n, r=1:K
            z_val[i,r] = callback_value(cb_data, z[i,r])
        end
        
        z_val = round.(z_val[:,:])
        z_val[map(v -> isnan(v) , z_val)] .= 0 # replace NaN values with 0
        z_val = Int.(z_val)

        # Calculate objective function and approximation difference
        pwl_obj = - Cnst + 0.5 * sum( x_val[i,j,r,s] for r=1:K, s=1:K, i=1:n, j=1:n)
        nl_obj = calculateObjective(dataset, w_val, z_val)

        # Check that the linearized objective value is less than or equal to the real objective value
        @assert nl_obj >= pwl_obj "Something weird inside the lazy constraint. PWL objective should be < NL obj. Got PWL obj: $pwl_obj and NL obj: $nl_obj"
        obj_relative_diff = (nl_obj - pwl_obj) / pwl_obj
        println("Lazy count $lazycount")
        println("  PWL obj: $(@sprintf("%.10f", pwl_obj))")
        println("  NL obj: $(@sprintf("%.10f", nl_obj))")
        println("  obj relative diff: $(@sprintf("%.10e", obj_relative_diff))") 
        println("  z_val : $z_val")
        println("  w_val : $w_val")

        # Allow for some impreciseness in the solution
        if abs(obj_relative_diff) >= TOL
            lazycount += 1;
            # Calculate new approximation
            w_cut = w_val[0 .< w_val .< Inf] # Only finite values
            size_w = length(w_cut)

            for i=1:n, j=1:n, p=1:size_w
                if (A[i,j] != 0)
                    a = -A[i,j]/w_cut[p] + ((k[i]*k[j])/(2*m)) # Slope
                    a = max(-1e8, min(1e8, a)) # Clip values for the slope
                    b = A[i,j]*(1 - log(w_cut[p])) # Intercept of each line
                    
                    # Submit constraints
                    for r=1:K, s=1:K
                        rhs = (M_upper[i,j,r,s] - b)
                        lhs = -x_val[i,j,r,s] + M_upper[i,j,r,s]*y_val[i,j,r,s] + a*w_val[r,s]
                        
                        if (!( lhs <= rhs )) 
                            lazycon = @build_constraint(a*w[r,s] + b - M_upper[i,j,r,s]*(1 - y[i,j,r,s]) <= x[i,j,r,s] ); # Lazy constraint
                            MOI.submit(model, MOI.LazyConstraint(cb_data), lazycon)
                        end
                    end
                end
            end
        end
    end # End of callback function
    
    # Tell JuMP/Gurobi/CPLEX to use our callback function
    MOI.set(model, MOI.LazyConstraintCallback(), lazy_callback_function)
    # ------------------------------------------------------------------------


    # --- Solve optimization model ------------------------------------------- 
    optimize!(model) # Solve 
    status      = string(termination_status(model))
    solvetime   = solve_time(model)
    obj_lb      = objective_bound(model);
    has_vals    = has_values(model)

    # Get objective value if there is a solution
    if has_vals
        obj_ub = objective_value(model);
    else
        obj_ub = nothing
    end
    
    # Get node count
    nodecount = Int64(MOI.get(model, MOI.NodeCount())) 
    # ------------------------------------------------------------------------


    # --- Recover variable values --------------------------------------------
    if has_vals
        w_opt = value.(w);
        x_opt = value.(x);
        z_opt = value.(z);
    else
        throw("No solution found. Revise the model output")
    end

    # Convert z to int 
    z_opt = round.(z_opt[:,:])
    z_opt[map(v -> isnan(v) , z_opt)] .= 0 # replace NaN values with 0
    z_opt = Int.(z_opt)

    # Check the tolerance in the objective value relative difference
    nonlinear_obj = calculateObjective(dataset, w_opt, z_opt)
    obj_relative_diff = (nonlinear_obj - obj_ub) / obj_ub
    if status == "OPTIMAL"
        @assert (obj_relative_diff <= TOL) "Tolerance of the objective value violated. Obj. value: $(@sprintf("%.10f", obj_ub)); NL obj: $(@sprintf("%.10f", nonlinear_obj)); Obj. relative diff: $(@sprintf("%.10e", obj_relative_diff))"
    end
    # ------------------------------------------------------------------------
    # SBM
    sbm = SBM(w_opt, "poisson") # assumes a poisson distribution

    # Print results
    println()
    println("w_opt =")
    display(w_opt)
    println()
    println("z_opt =")
    display(z_opt)
    println()
    println("Solve time:                $solvetime")
    println("Status:                    $status")
    println("Objective value:           $obj_ub")
    println("Objective lower bound:     $obj_lb")
    println("Non-linear objective:      $(@sprintf("%.10f", nonlinear_obj))")
    println("Obj relative diff:         $(@sprintf("%.10e", obj_relative_diff))")
    println("Nodes processed:           $nodecount")
    println("Lazy constraints:          $lazycount")

    # Build OptResults
    iterations = nothing
    opt_results = OptResults(obj_lb, obj_ub, status, solvetime, iterations, nodecount, lazycount, num_cons)
    displayResults(opt_results)
    return sbm, z_opt, opt_results
end



function to_standard_form(x::Matrix{Int})::Matrix{Int}
    n,K = size(x)
    
    mapx = zeros(Int, K)
    new_x = zeros(Int, n, K)

    gmax = 0

    for i=1:n
        r = findfirst(x[i,:] .== 1)
        if mapx[r] == 0
            gmax += 1
            mapx[r] = gmax
        end
        new_x[i,mapx[r]] = 1
    end

    return new_x
end

function randomAssignments(dataset::Dataset; seed::Union{Nothing, Int}=nothing)::Matrix{Int}
    """ Generates random assignments of n nodes to K groups

    Parameters
    ----------
    dataset : Dataset
        Dataset representing an observed graph.

    Returns
    -------
    z : Array{Int64,2}
        Matrix of assignments of nodes to groups
    """
    n = dataset.n
    K = dataset.n_communities
    if isnothing(seed)
        Random.seed!(seed)
    end

    assign = rand(1:K, n)

    z = zeros(Int, n, K)
    for i=1:n
        z[i,assign[i]] = 1
    end
    @assert sum(z) == n
    return z
end

function maximizationStep(dataset::Dataset, z::Matrix{Int})::Matrix{Float64}
    """ Calculates the optimal probability matrix given the assignments

    Parameters
    ----------
    dataset : Dataset
        Dataset representing an observed graph.
    z : Array{Int64,2}
        Matrix of assignments of nodes to groups

    Returns
    -------
    w : Array{Int64,2}
        Matrix of probabilities
    """
    A = dataset.A
    n = dataset.n
    m = dataset.m
    K = dataset.n_communities
    k = dataset.k
    if ~((n,K) == size(z))
        throw(ArgumentError("The dimensions of the assignment matrix do not match the dimensions of the dataset."))
    end

    w = zeros(K,K)
    for g=1:K, h=1:K
        numerator = 0.5*sum(A[i,j]*z[i,g]*z[j,h] for i=1:n, j=1:n)
        denominator = 0.5*sum( ((k[i]*k[j])/(2*m)) * z[i,g]*z[j,h] for i=1:n, j=1:n)
        if numerator == 0
            w[g,h] = 0
        else
            w[g,h] = numerator / denominator
        end
    end
    return w
end



function expectationStep(dataset::Dataset, 
                         w::Matrix{Float64}, 
                         z::Matrix{Int},
                         variant::String,
                         time_limit::Float64)::Matrix{Int}
    if variant in ["ls1","ls2"]
        return expectationStep_LS(dataset, w, z, variant, time_limit)
    elseif variant == "exact"
        return expectationStep_exact(dataset, w, z, time_limit)
    end
end


function expectationStep_LS(dataset::Dataset, 
                            w::Matrix{Float64}, 
                            z::Matrix{Int},
                            variant::String, 
                            time_limit::Float64)::Matrix{Int}

    # Check input argument 'e_step'
    if ~(variant in ["ls1","ls2"])
        throw(ArgumentError(string("Invalid variant for E-step: ", variant)))
    end

    n = dataset.n
    K = dataset.n_communities

    start = time(); # Start counting time

    objval = calculateObjective(dataset, w, z)

    improved = true
    while true

        ## Base cases
        # Break if there is no improvement
        if ~improved
            break # Local optimum
        end
        improved = false

        for i in Random.randperm(n)
            for g in Random.randperm(K)
                # Break if time limit is exceeded
                availableTime = time_limit - (time() - start)
                if availableTime <= 0
                    return z
                end
                
                # Consider z' constructed from z by relocating vertex i to community g
                new_z = copy(z) 
                new_z[i,:] = zeros(Int, K)
                new_z[i,g] = 1

                new_w = copy(w)

                # Integrated maximization step
                if variant == "ls2"
                    new_w = maximizationStep(dataset, new_z)
                end

                # Re-calculate objective value
                new_objval = calculateObjective(dataset, new_w, new_z)

                if new_objval < objval
                    w = new_w
                    z = new_z
                    objval = new_objval

                    improved = true
                end
            end
        end
    end

    return z
end



function expectationStep_exact(dataset::Dataset, 
                               w::Matrix{Float64}, 
                               z::Matrix{Int},
                               time_limit::Float64,
                               TOL::Float64=1e-6)::Matrix{Int}

    A = dataset.A
    K = dataset.n_communities
    n = dataset.n
    m = dataset.m
    k = dataset.k

    verbose = true
    symm_break_const = true

    w[w .== 0] .= 1e-12
    W = zeros(n,n,K,K);
    for i=1:n, j=1:n, g=1:K, h=1:K
        if A[i,j] != 0
            W[i,j,g,h] = 0.5*( -A[i,j]*log(w[g,h]) + ((k[i]*k[j])/(2*m))*w[g,h] )
        else
            W[i,j,g,h] = 0.5*((k[i]*k[j])/(2*m))*w[g,h]
        end
    end

    solver = CPLEX.Optimizer
    model = Model(solver)
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", verbose)
    set_optimizer_attribute(model, "CPX_PARAM_TILIM", time_limit)
    set_optimizer_attribute(model, "CPX_PARAM_EPGAP", TOL)
    set_optimizer_attribute(model, "CPX_PARAM_EPINT", 1e-9) # integrality tolerance
    set_optimizer_attribute(model, "CPX_PARAM_EPRHS", 1e-9) # feasibility tolerance
    set_optimizer_attribute(model, "CPX_PARAM_EPOPT", 1e-9) # optimality tolerance
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", false) 
    
    # Variables
    @variable(model, z[1:n,1:K], Bin);

    # Objective
    @objective(model, Min, sum(
            (W[i,j,g,h] * z[i,g] * z[j,h])
            for g=1:K, h=1:K, i=1:n, j=1:n));

    # Assignment constraints
    @constraint(model, assign[i = 1:n], sum(z[i,g] for g=1:K) == 1); # assignment constraint
    
    # Symmetry breaking constraints
    if symm_break_const == true
        @constraint(model, conObj1, z[1,1] == 1); # object 1 must be in cluster 1
        @constraint(model, symBreak[r=2:(K-1), j=r:n], sum(z[i,l] for l=1:(r-1), i=2:(j-1)) - sum(z[j,l] for l=1:r) <= j - 3 )
    end

    # Solve
    start = time()
    optimize!(model) 
    status = string(termination_status(model))
    println("status = $status")
    solvetime = time() - start

    obj_lb = objective_bound(model);
    if isnan(obj_lb)
        obj_lb = -Inf
    end
    obj_ub = objective_value(model);
    if isnan(obj_ub)
        obj_ub = Inf
        return z
    end

    ############################
    # Recover variable values
    z = value.(z)
    z = round.(z[:,:])
    z[map(v -> isnan(v) , z)] .= 0 # replace NaN values with 0
    z = Int.(z)

    return z
end



function expectationMaximization(dataset::Dataset, 
                                 e_step::String;
                                 time_limit::Float64=600.0,
                                 seed::Union{Nothing, Int}=nothing,
                                 z::Union{Nothing,Matrix{Int}}=nothing)::Tuple{SBM, Matrix{Int}, OptResults}

    # Check input argument 'e_step'
    if ~(e_step in ["ls1","ls2","exact"])
        throw(ArgumentError(string("Invalid argument for E-step: ", e_step)))
    end

    # Start counting time
    start = time(); 

    # Initialize with random assignments
    if isnothing(z)
        z = randomAssignments(dataset, seed=seed)
    end
    w = maximizationStep(dataset, z)

    # Calculate objective value
    objval = calculateObjective(dataset, w, z)

    status = nothing
    iterations = 0
    improved = true
    while true
        ## Base cases
        # Break if there is no improvement
        if ~improved
            status = "LocalOptimum"
            break
        end
        improved = false

        # Break if time limit is exceeded
        availableTime = time_limit - (time() - start)
        if availableTime <= 0
            status = "TIME_LIMIT"
            break
        end

        # Expectation Step
        z = expectationStep(dataset, w, z, e_step, availableTime)
        z = to_standard_form(z)

        # Maximization Step
        w = maximizationStep(dataset, z)

        # Check if there was an improvement (otherwise it has reached a local optimum)
        new_objval = calculateObjective(dataset, w, z)
        improved = (new_objval < objval) # Check for improvement
        objval = new_objval # Update objective value (it does not deteriorate in the EM)
        
        # Increment iterations
        iterations += 1
    end
    solvetime = time() - start;

    # SBM
    sbm = SBM(w, "poisson") # assumes a poisson distribution

    # Build OptResults
    obj_lb = -Inf
    obj_ub = objval
    nodecount = nothing
    lazycount = nothing
    num_constraints = nothing
    opt_results = OptResults(obj_lb, obj_ub, status, solvetime, iterations, nodecount, lazycount, num_constraints)
    displayResults(opt_results)
    return sbm, z, opt_results
end



function multiStartEM(dataset::Dataset, 
                      e_step::String;
                      time_limit::Float64=600.0,
                      n_trials::Union{Nothing,Int}=nothing)

    # Check input argument 'e_step'
    if ~(e_step in ["ls1","ls2","exact"])
        throw(ArgumentError(string("Invalid argument for E-step: ", e_step)))
    end

    # Start counting time
    start = time(); 

    best_ub = Inf
    trial = 1
    sbm = nothing
    z = nothing
    status = nothing
    while true
        # Maximum number of trials
        if (~isnothing(n_trials) && trial > n_trials)
            status = "LocalOptimum"
            break
        end

        # Break if time limit is exceeded
        availableTime = time_limit - (time() - start)
        if availableTime <= 0
            status = "TIME_LIMIT"
            break
        end

        sbm_t, z_t, res_t = expectationMaximization(dataset, e_step, time_limit)

        # Check if there was an improvement
        if res_t.UB < best_ub
            best_ub = res_t.UB # Update best solution
            sbm = sbm_t
            z = z_t
        end

        trial += 1
    end

    # Build OptResults
    obj_lb = -Inf
    obj_ub = best_ub
    solvetime = (time() - start)
    nodecount = nothing
    lazycount = nothing
    num_constraints = nothing
    opt_results = OptResults(obj_lb, obj_ub, status, solvetime, trial, nodecount, lazycount, num_constraints)
    displayResults(opt_results)
    return sbm, z, opt_results
end

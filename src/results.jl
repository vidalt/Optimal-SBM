struct OptResults
    LB::Float64                         # Lower Bound on the objective
    UB::Float64                         # Upper Bound on the objective
    status::String                      # Status
    solvetime::Float64                  # Solve time in seconds
    iterations::Union{Nothing,Int}      # Number of iterations (for local search)
    nodecount::Union{Nothing,Int}       # Node count of the Branch-and-Bound tree
    lazycount::Union{Nothing,Int}       # Number of lazy constraints added during solving
    num_constraints::Union{Nothing,Int} # Number of constraints in the original model (before solving)

    function OptResults(lb::Float64, ub::Float64, status::String, solvetime::Float64,
                        iterations::Union{Nothing,Int}, nodecount::Union{Nothing,Int}, 
                        lazycount::Union{Nothing,Int}, num_constraints::Union{Nothing,Int})
        if ~(lb <= ub) && ~(isapprox(lb, ub, rtol=1e-4))
            throw(ArgumentError("Lower Bound value must be less than or equal to Upper Bound. Got values LB=$lb and UB=$ub"))
        end
        if (status == "OPTIMAL") && !(isapprox(lb, ub, rtol=1e-4))
            throw(ArgumentError("An optimal solution implies that the Lower Bound is equal to the Upper Bound. Got values LB=$lb and UB=$ub"))
        end
        if solvetime < 0
            throw(DomainError("Solve time must be a non-negative number (in seconds)."))
        end
        if (!isnothing(iterations)) && (iterations < 0)
            throw(DomainError("Number of iterations must be a positive integer number."))
        end
        if (!isnothing(nodecount)) && (nodecount < 0)
            throw(DomainError("Number of nodes in B&B must be a non-negative integer number."))
        end
        if (!isnothing(lazycount)) && (lazycount < 0)
            throw(DomainError("Number of lazy constraints added must be a non-negative integer number."))
        end
        if (!isnothing(num_constraints)) && (num_constraints < 0)
            throw(DomainError("Number of constraints must be a non-negative integer number."))
        end

        return new(lb, ub, status, solvetime, iterations, nodecount, lazycount, num_constraints)
    end
end

function displayResults(opt_results::OptResults)
    """
    Method to display optimization results.

    Parameters
    ----------
    opt_results : OptResults
        Optimization results object.
    """
    print_string = string("\n",
          "--------------- Opt Results ---------------\n",
          "Obj. LB = ",round(opt_results.LB, digits=5),"\n",
          "Obj. UB = ",round(opt_results.UB, digits=5),"\n",
          "Status: $(opt_results.status)\n",
          "Solve time: ",round(opt_results.solvetime, digits=5)," seconds\n")
    if !isnothing(opt_results.iterations)
        print_string = string(print_string, "Iterations: $(opt_results.iterations)\n")
    end
    if !isnothing(opt_results.nodecount)
        print_string = string(print_string, "Number of nodes: $(opt_results.nodecount)\n")
    end
    if !isnothing(opt_results.num_constraints)
        print_string = string(print_string, "Number of constraints: $(opt_results.num_constraints)\n")
    end
    if !isnothing(opt_results.lazycount)
        print_string = string(print_string, "Number of lazy constraints: $(opt_results.lazycount)\n")
    end
    print(print_string)
end

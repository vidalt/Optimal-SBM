struct SBM
    w::Matrix{Float64}      # Matrix of probabilities
    dist::String            # Distribution (either "poisson" or "bernoulli")
    n_communities::Int      # Number of communities/clusters
    variant::String

    # Constructor function for struct SBM
    function SBM(w::Matrix{Float64}, dist::String, variant::String)
        # Check if w is a square matrix
        if size(w,1) != size(w,2)
            throw(ArgumentError("Matrix w must be a square matrix."))
        end
        # Check if elements in w are strictly positive
        if any(w .< 0)
            throw(DomainError("All elements in matrix of probabilities w must be strictly positive. Got w = $w"))
        end
        # Check if matrix is symmetric
        if ~isapprox(w, w', rtol=1e-4) # TODO: define global variable for rtol
            throw(ArgumentError("Matrix w must be symmetric. Got w = $w"))
        end
        q = size(w,1)

        # Check input argument 'distribution'
        if ~(dist in ["poisson","bernoulli"])
            throw(ArgumentError(string("Invalid distribution: ", dist)))
        end
        if dist == "bernoulli"
            # Check if w is a matrix of 0-1 values
            if any(w .> 1)
                throw(DomainError("When using a bernoulli distribution all elements in matrix of probabilities w must be between 0 and 1. Got w = $w"))
            end
        end
        if ~(variant in ["standard","degree-corrected"])
            throw(ArgumentError("Invalid SBM variant: $variant"))
        end
        return new(w, dist, q, variant)
    end

    function SBM(w::Matrix{Float64}, dist::String)
        variant = "standard"
        return SBM(w, dist, variant)
    end
end

function generate(sbm::SBM, n_per_community::Vector{Int}, seed::Union{Nothing,Int}=nothing)::Dataset
    """
    Generates a graph from the Stochastic Block Model.

    Parameters
    ----------
    sbm : SBM
        An SBM structure.
    n_per_community : Vector{Int}
        Number of nodes in each community.
    seed : Union{Nothing,Int}
        Seed for the random number generator.

    Returns
    -------
    adj_matrix : Matrix{Int}
        Adjacency matrix of the generated graph
    """
    if sbm.variant != "standard"
        throw("not implemented")
    end

    # Check input argument 'n_per_community'
    if any(n_per_community .<= 0)
        throw(DomainError(string("Invalid value in argument n_per_community: ",
                                 "the number of nodes in each community must be strictly positive.")))
    end
    # Check if n_per_community matches the dimensions of w
    if ~(length(n_per_community) == sbm.n_communities)
        throw(ArgumentError("Number of communities in n_per_community must match the number of communities in SBM"))
    end

    # Set the random seed if one was provided
    if seed != nothing
        Random.seed!(seed)
    end
    # Define number of vertices, initialize adjacency matrix and community labels
    n_vertices = sum(n_per_community);
    adj_matrix = zeros(Int, n_vertices,n_vertices);
    community_labels = convert(Array{Int}, vcat([i*ones(ni) for (i,ni) in enumerate(n_per_community)]...));
    # Construct adjacency matrix
    for i=1:n_vertices
        for j=1:i
            proba = sbm.w[community_labels[i], community_labels[j]];
            if sbm.dist == "poisson"
                if i == j
                    # Note that this implies that the expected number of self-edges at a vertex in group r is 1/2 * ωrr 
                    # because of the factor of two in the definition of the diagonal elements of the adjacency matrix.
                    proba = proba / 2 # "The factor of half is included solely because it makes the algebra easier"
                end
                pois = Distributions.Poisson(proba);
                adj_matrix[i,j] = rand(pois);
            elseif sbm.dist == "bernoulli"
                if i == j
                    adj_matrix[i,j] = 0;
                else
                    bernoulli = Bernoulli(proba);
                    adj_matrix[i,j] = rand(bernoulli);
                end
            end
        end
    end
    # "We have adopted the common convention that a self-edge is represented by A[i,j] = 2
    # (and not 1 as one might at first imagine)"
    adj_matrix += transpose(adj_matrix) #- Diagonal(adj_matrix)

    # Instantiate a dataset
    n_communities = length(n_per_community)
    dataset = Dataset(adj_matrix, n_communities)
    return dataset
end

function generate(probability_matrix::Matrix{Float64}, n_per_community::Vector{Int};
                            distribution::String="poisson", seed::Union{Nothing, Int}=nothing)::Dataset
    """
    Generates a graph from the Stochastic Block Model.

    Parameters
    ----------
    probability_matrix : Matrix{Float64}
        Matrix of expected number of edges under Poisson distribution.
    n_per_community : Vector{Int}
        Number of nodes in each community.
    distribution : String
        Distribution to be used for generating edges (either "poisson" or "bernoulli").
    seed : Union{Nothing,Int}
        Seed for the random number generator.

    Returns
    -------
    dataset : Dataset
        Dataset representing a graph.
    """
    sbm = SBM(probability_matrix, distribution)
    return generate(sbm, n_per_community, seed)
end




function getObjectiveConstant(dataset::Dataset)::Float64
    """ Calculate constant term of the objective function """
    A = dataset.A
    n = dataset.n
    m = dataset.m
    k = dataset.k

    Cnst = 0
    for i=1:n, j=1:(i-1)
        if A[i,j] == 0
            Cnst += -1
        else
            Cnst += ( A[i,j] * log(k[i]*k[j] / (2*m)) - log(factorial(A[i,j])) )
        end
    end

    for i=1:n
        if A[i,i] == 0
            Cnst += -1
        else
            Cnst += ( 0.5 * A[i,i] * log((0.5*k[i]^2)/(2*m)) - log(0.5*factorial(A[i,i])))
        end
    end
    return Cnst
end


function getWUpperBound(dataset::Dataset)::Float64
    # Upper bound on w
    n = dataset.n
    m = dataset.m
    A = dataset.A
    k = dataset.k

    rho = 0
    for i=1:n, j=1:n
        if A[i,j] != 0
            rho = maximum([rho, (A[i,j] / (k[i]*k[j]))])
        end
    end
    w_upper = 2 * m * rho
    return w_upper
end


function calculateObjective(dataset::Dataset, w::Matrix{Float64}, z::Matrix{Int})::Float64
    """ Calculates the maximum log-likelihood value

    Parameters
    ----------
    dataset : Dataset
        Dataset representing an observed graph.
    sbm : SBM
        An Stochastic Block Model
    z : Matrix{Int}
        Matrix of assignments of nodes to groups

    Returns
    -------
    L : Float64
        Maximum log-likelihood value
    """
    A = dataset.A
    n = dataset.n
    m = dataset.m
    q = dataset.n_communities
    k = dataset.k
    if ~((n,q) == size(z))
        throw(ArgumentError("The dimensions of the assignment matrix do not match the dimensions of the dataset."))
    end

    # Get objective constant
    Cnst = getObjectiveConstant(dataset)

    # Calculate objective value
    L = 0
    for g=1:q, h=1:q, i=1:n, j=1:n
        if z[i,g] * z[j,h] == 1
            if A[i,j] != 0
                L += 0.5*(- A[i,j]*log(w[g,h]) + ((k[i]*k[j])/(2*m)) * w[g,h])
            else
                L += 0.5*( ((k[i]*k[j])/(2*m)) * w[g,h] )
            end
        end
    end
    L = L - Cnst
    return L
end


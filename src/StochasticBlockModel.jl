module StochasticBlockModel

using Distributions, Random, JuMP, AmplNLWriter, CPLEX, Printf

include("datasets.jl")
include("sbm.jl")
include("results.jl")
include("exact_methods.jl")
include("minlp_ampl.jl")
include("expectation_maximization.jl")

export generate, Dataset

end # module

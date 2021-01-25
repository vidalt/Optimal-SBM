using DataFrames
using CSV
using DelimitedFiles
using Statistics
using Glob
using LinearAlgebra
include("../src/StochasticBlockModel.jl")

all_datasets = glob("instances/S1/in/*.in")
output_path = "experiments/experiment-S1-heuristics.txt"

heuristic_methods = ["ls1","ls2","exact"]

time_limit = 20.0
N_trials = 50

df = DataFrame(instance_path = String[],
               dataset_name = String[],
               n = Int[],
               m = Int[],
               q = Int[],
               avg_degree = Float64[],  
               method = String[],
               LB = Float64[],
               UB = Float64[],
               status = String[],
               solvetime = Float64[],
               trial = Int[], 
               seed = Int[],
               strong_assort = Bool[],
               weak_assort = Bool[],
               z = Matrix{Int}[],
               w = Matrix{Float64}[] 
)

# Write column names of DataFrame to output file
df_columns = string(join(string.(names(df)), "\t"),"\n")
open(output_path, "a+") do io
   write(io, df_columns)
end

for filepath in all_datasets
    dataset_path = split(filepath, "/")[end]
    println(filepath)

    # Load Dataset
    dataset = StochasticBlockModel.Dataset(filepath)

    for method in heuristic_methods
        println("method: $method")
        for trial in 1:N_trials
            println("Trial: $trial")
            # Sample a random seed
            seed = abs(rand(Int))
            println("seed: $seed")
            
            # Estimate SBM
            sbm, z, opt_results = StochasticBlockModel.expectationMaximization(dataset, 
                                                                               method, 
                                                                               time_limit=time_limit, 
                                                                               seed=seed, 
                                                                               z=nothing)
            
            # Average degree
            avg_degree = mean(dataset.k)

            # Strong assortativity
            strong_assort = (minimum(diag(sbm.w)) >= maximum([sbm.w[i,j] 
                             for i=1:dataset.n_communities, j=1:dataset.n_communities if j>i]))

            # Weak assortativity
            weak_assort = all([sbm.w[i,i] >= maximum([sbm.w[i,j] for j=1:dataset.n_communities if i!=j]) 
                               for i=1:dataset.n_communities])

            # Append to data frame
            row = (filepath, dataset_path, dataset.n, dataset.m, dataset.n_communities, avg_degree, method,
                   opt_results.LB, opt_results.UB, string(opt_results.status), opt_results.solvetime,
                   trial, seed, strong_assort, weak_assort, z, sbm.w)
            push!(df, row)

            # Write new line to output file
            CSV.write(output_path, df[end:end,:], delim="\t", append=true)

            # Save prob matrix w and assignments z
            res_folder = joinpath(join(split(output_path,"/")[1:end-1], "/"), string("n",dataset.n))
            if ~isdir(res_folder)
                mkdir(res_folder)
            end
            
            # Save w
            w_savepath = replace(joinpath(res_folder, dataset_path),".in"=>"-t$trial-$method.w")
            writedlm(w_savepath, sbm.w, '\t')

            # Save z
            z_savepath = replace(joinpath(res_folder, dataset_path),".in"=>"-t$trial-$method.z")
            writedlm(z_savepath, z, '\t')
        end
    end
end
using DataFrames
using CSV
using DelimitedFiles
using Glob
using Statistics
using LinearAlgebra
include("../src/StochasticBlockModel.jl")

println("Start execution")

println(pwd())
time_limit = 20.0
output_path = "experiments/experiment-S1-milp.txt"
verbose = false
file_num_start = 0

### -------------------------------------------------------------------
all_datasets = glob("instances/S1/in/*.in")

println(string("Datasets found: ", length(all_datasets)))

symm_break_const_list = [true]
symm_w_const_list = [true]

### -------------------------------------------------------------------

method = "milp"

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
               nodecount = Int[],
               lazycount = Int[],
               symm_break_const = Bool[],
               symm_w_const = Bool[],
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

for (num_file, filepath) in enumerate(all_datasets)
    if num_file < file_num_start
        continue
    end

    dataset_path = split(filepath, "/")[end]
    println("File number: $num_file")
    println(filepath)

    # Load Dataset
    dataset = StochasticBlockModel.Dataset(filepath)
    println(string("Size: ", dataset.n))


    # enforce symmetry on w
    for symm_w in symm_w_const_list
        # symmetry breaking constraints
        for symm_break in symm_break_const_list
            # Estimate SBM
            sbm, z, opt_results = StochasticBlockModel.MILP(dataset, time_limit=time_limit, symm_break_const=symm_break, symm_w_const=symm_w)
            
            avg_degree = mean(dataset.k)
            nodecount = opt_results.nodecount
            if isnothing(nodecount)
                nodecount = 0
            end
            lazycount = opt_results.lazycount
            if isnothing(lazycount)
                lazycount = 0
            end

            # Strong assortativity
            strong_assort = (minimum(diag(sbm.w)) >= maximum([sbm.w[i,j] 
                            for i=1:dataset.n_communities, j=1:dataset.n_communities if j>i]))

            # Weak assortativity
            weak_assort = all([sbm.w[i,i] >= maximum([sbm.w[i,j] for j=1:dataset.n_communities if i!=j]) 
                                                    for i=1:dataset.n_communities])

            row = (filepath, dataset_path, dataset.n, dataset.m, dataset.n_communities, avg_degree, method,
                    opt_results.LB, opt_results.UB, string(opt_results.status), opt_results.solvetime,
                    nodecount, lazycount, symm_break, symm_w, strong_assort, 
                    weak_assort, z, sbm.w)

            push!(df, row)

            # Write new line to output file
            CSV.write(output_path, df[end:end,:], delim="\t", append=true)

            # Save prob matrix w and assignments x
            res_folder = joinpath(join(split(output_path,"/")[1:end-1], "/"), string("n",dataset.n))
            if ~isdir(res_folder)
                mkdir(res_folder)
            end
            
            # Save w
            w_savepath = replace(joinpath(res_folder, dataset_path),".in"=>"-$method.w")
            writedlm(w_savepath, sbm.w, '\t')

            # Save z
            z_savepath = replace(joinpath(res_folder, dataset_path),".in"=>"-$method.z")
            writedlm(z_savepath, z, '\t')
        end
    end
end

println("End of execution.")

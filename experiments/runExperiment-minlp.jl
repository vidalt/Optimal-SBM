using Glob
using Statistics
using LinearAlgebra
include("../src/StochasticBlockModel.jl")

println("Start execution")

println(pwd())
file_num_start = 0

### -------------------------------------------------------------------
all_datasets = glob("instances/S1/in/n8-*.in")[1:10]

out_folder = "experiments/"

ampl_path   = "/opt/ampl.linux64/ampl"

symm_break_const_list = [true]
symm_w_const_list = [true]

### -------------------------------------------------------------------

method = "minlp"

mod_file = "src/minlp.mod"


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

    dataset_name = split(basename(filepath),".")[1]

    StochasticBlockModel.MINLP_AMPL(dataset, dataset_name, mod_file, out_folder, ampl_path)

end

println("End of execution.")

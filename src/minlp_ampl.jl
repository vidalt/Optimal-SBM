

### Create data file
function write_data(output_path, dataset, Cnst, w_upper)
    
    open(output_path, "w") do io
        write(io, "param n := $(dataset.n);\n")
        write(io, "param m := $(dataset.m);\n")
        write(io, "param K := $(dataset.n_communities);\n")
        write(io, "\n")
        write(io, "param\tA:")

        for i=1:dataset.n
            write(io, "\t$i")
        end
        write(io, ":=")
        for i=1:dataset.n
            write(io, "\n\t\t$i")
            for j=1:dataset.n
                write(io, "\t$(dataset.A[i,j])")
            end
        end
        write(io, ";\n")

        write(io, "param\tk :=")
        for i=1:dataset.n
            write(io, "\n\t\t$i\t$(dataset.k[i])")
        end
        write(io, ";\n")
        
        write(io, "param\tCnst := $Cnst;\n")
        write(io, "param\tw_upper := $w_upper;\n")
    end;
end


### Create run file
function write_run_file(run_file, mod_file, data_file, log_file)
    open(run_file, "w+") do io
        write(io, "reset;\n")
        write(io, "model $mod_file;\n")
        write(io, "data $data_file;\n")
        write(io, "option log_file '$log_file';\n")
        write(io, "expand L, Constraint;\n")
        write(io, "option solver couenne;\n")
        write(io, "solve;\n")
        write(io, "display n;\n")
        write(io, "display m;\n")
        write(io, "display K;\n")
        write(io, "display A;\n")
        write(io, "display k;\n")
        write(io, "display Cnst;\n")
        write(io, "display w;\n")
        write(io, "display z;\n")
        write(io, "display L;\n")
    end
end



### Read log file and extract information
function extract_solution(dataset_name, log_file, sol_file)
    io = open(log_file)

    n = m = K = -1
    lincuts_root = lincuts_total = separation_time = total_solve_time = bb_time =\
        obj_lb = obj_ub = gap = bb_nodes = cnst = -1
    z = w = nothing
    z_mat = nothing  # Variable to track if z is written in a matrix or vector format
    status = "?"
    processing_x = false
    processing_w = false

    while !eof(io)
        line = readline(io)

        if startswith(line, "n = ")
            n = parse(Int, strip(split(line, "n = ")[2]))
        end
        if startswith(line, "m = ")
            m = parse(Int, strip(split(line, "m = ")[2]))
        end 
        if startswith(line, "K = ")
            K = parse(Int, strip(split(line, "K = ")[2]))
        end 
        if startswith(line, "couenne: Optimal")
            status = "OPTIMAL"
        end
        if startswith(line, "couenne: Optimization interrupted on limit.")
            status = "TIME_LIMIT"
        end
        if startswith(line, "Linearization cuts added at root node:")
            lincuts_root = parse(Int, strip(split(line, ":")[2]))
        end
        if startswith(line, "Linearization cuts added in total:")
            lincuts_total = parse(Int, strip(split(split(line, ":")[2], "(")[1]))
        end 
        if startswith(line, "Linearization cuts added in total:")
            separation_time = parse(Float64, strip(split(split(line, "separation time:")[2], "s)")[1]))
        end 
        if startswith(line, "Total solve time:")
            total_solve_time = parse(Float64, strip(split(split(line, ":")[2], "(")[1])[1:end-1])
        end
        if startswith(line, "Total solve time:")
            bb_time = parse(Float64, strip(split(split(line, "(")[2], "s in branch-and-bound)")[1]))
        end 
        if startswith(line, "Lower bound:")
            obj_lb = parse(Float64, strip(split(line, ":")[2]))
        end
        if startswith(line, "Upper bound:")
            obj_ub = parse(Float64, strip(split(split(line, ":")[2], "(")[1]))
        end 
        if startswith(line, "Upper bound:")
            gap = strip(split(split(line, "(gap: ")[2], ")")[1])
            if gap == "--"
                gap = NaN
            else
                gap = parse(Float64, gap[1:end-1])
            end
        end
        if startswith(line, "Branch-and-bound nodes:")
            bb_nodes = parse(Int, strip(split(line, ":")[2]))
        end
        if startswith(line, ";") && processing_x
            processing_x = false
        end
        if processing_x
            if z_mat
                if startswith(line, ":")
                else
                    elements = parse.(Float64, split(strip(line)))
                    println(elements)
                    i = Int(elements[1])
                    z[i,:] = elements[2:end]
                end
            else
                i,j,val = parse.(Float64, split(strip(line)))
                i = Int(i)
                j = Int(j)
                z[i,j] = val
            end
        end
        
        if startswith(line, "z :=")
            z = zeros(Float64, n, K)
            processing_x = true
            z_mat = false
        end 
        if startswith(line, "z [*,*]")
            z = zeros(Float64, n, K)
            processing_x = true
            z_mat = true
        end 
        if startswith(line, ";") && processing_w
            processing_w = false        
        end
        if processing_w
            i,g,val = parse.(Float64, split(strip(line)))
            i = Int(i)
            g = Int(g)
            w[i,g] = val
        end
        if startswith(line, "w :=")
            w = zeros(Float64, K, K)
            processing_w = true
        end 
        if startswith(line, "Cnst =")
            cnst = parse(Float64, strip(split(line, "=")[2]))
        end
    end
    close(io)
    
    println("Status: $status")
    println("Linearization cuts at root: $lincuts_root")
    println("Linearization cuts total: $lincuts_total")
    println("Separation time: $separation_time")
    println("Total solve time: $total_solve_time")
    println("B&B time: $bb_time")
    println("LB: $obj_lb")
    println("UB: $obj_ub")
    println("gap: $gap")
    println("B&B nodes: $bb_nodes")
    println("z = ")
    println(z)
    println("w = ")
    println(w)
    println("Cnst = $cnst")
 
    # Convert array type to Int
    z_int = convert.(Int, z)

    # Save solution to file
    open(sol_file, "w") do io
        write(io, "\tSolution\n")
        write(io, "dataset_name\t$dataset_name\n")
        write(io, "n\t$n\n")
        write(io, "m\t$m\n")
        write(io, "K\t$K\n")
        write(io, "Cnst\t$cnst\n")
        write(io, "method\tminlp\n")
        write(io, "LB\t$obj_lb\n")
        write(io, "UB\t$obj_ub\n")
        write(io, "status\t$status\n")
        write(io, "solvetime\t$total_solve_time\n")
        write(io, "bb_solvetime\t$bb_time\n")
        write(io, "nodecount\t$bb_nodes\n")
        write(io, "lincuts_root\t$lincuts_root\n")
        write(io, "lincuts_total\t$lincuts_total\n")
        write(io, "separation_time\t$separation_time\n")
        write(io, "gap\t$gap\n")
        write(io, "w\t$w\n")
        write(io, "z\t$z_int\n")
        write(io, "z_float\t$z\n") # TODO: convert to int/float?
    end
end


function MINLP_AMPL(dataset, dataset_name, mod_file, out_folder, ampl_path)

    println("Processing dataset:", dataset_name)
    
    data_file = string(out_folder, dataset_name, ".dat")

    Cnst = getObjectiveConstant(dataset)
    w_upper = getWUpperBound(dataset)
    write_data(data_file, dataset, Cnst, w_upper)
    
    method = split(basename(mod_file), ".")[1]
    
    ### Create run file
    run_file = string(out_folder, method, "-", dataset_name, ".run")
    log_file = string(out_folder, method, "-", dataset_name, ".log")
    
    write_run_file(run_file, mod_file, data_file, log_file)
    
    ### Run ampl couenne
    println("Running couenne")
    run(`$ampl_path $run_file`)

    log_file = string(out_folder, method, "-", dataset_name, ".log")

    # Extract solution details from log
    sol_file = string(out_folder, method, "-", dataset_name, ".out")
    extract_solution(dataset_name, log_file, sol_file)
    println()
end


#!/usr/bin/env julia
# Python and Advanced Data Structures
# https://teaching.healthtech.dtu.dk/22110/index.php/Python_and_Advanced_Data_Structures

# Ex 1
function read_experiment!(outdict, f)
    for line = eachline(f)
        line_split = split(line, "\t")
        key = string(line_split[1])
        nums = [parse(Float64, x) for x in line_split[2:end]]
        outdict[key] = append!(get(outdict, key, Float64[]), nums)
    end
end

function read_experiments(out, files...)
    expdict = Dict{String,Array{Float64,1}}()
    for file in files
        open(f->read_experiment!(expdict, f), file)
    end
    for (key, value) in sort(collect(expdict))
        print(out, key, "\t")
        println(out, join(value, "\t"))
    end
end

open("data/test.dat", "w") do f
    read_experiments(f, "data/test1.dat", "data/test2.dat", "data/test3.dat")
end

# Ex 2

matrix = readdlm("data/matrix.dat", '\t', Int, '\n')
#matrix1 = readdlm("data/matrix1", '\t', Int, '\n')
#matrix2 = readdlm("data/matrix2", '\t', Int, '\n')
#dot = matrix1*matrix2
transpose(matrix)

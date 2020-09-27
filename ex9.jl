#!/usr/bin/env julia
# Comprehension and Generators
# https://teaching.healthtech.dtu.dk/22110/index.php/Comprehension_and_Generators

# Ex 1
matrix1 = readdlm("data/matrix1", '\t', Int, '\n')
matrix2 = readdlm("data/matrix2", '\t', Int, '\n')

dot = matrix1*matrix2
matrix_product(a,b) = [sum([a[i,r]*b[r,j] for r = 1:size(a,2)])
    for i = 1:size(a,1), j = 1:size(b,2)]

matrix_product(matrix1,matrix2)

# Ex 2
function _read_experiment(f, acc_match)
    status = Float64[]
    for line = eachline(f)
        line_split = split(line, "\t")
        acc = string(line_split[2])
        if line_split[1] == "COL_CLASSES"
            status = [parse(Int64, x) for x in line_split[4:end]]
        elseif acc == acc_match
            nums = [parse(Float64, x) for x in line_split[4:end]]
            return nums[status .== 0], nums[status .== 1]
        end
    end
    return Float64[], Float64[]
end


function read_experiment(acc_search, outfilename)
    out = open(outfilename, "w")
    cancer, control = open(f->_read_experiment(f, acc_search), "data/dna-array.dat")
    if length(cancer) == 0 && length(control) == 0
        println("Accession number not found.")
        close(out)
        return nothing
    end
    min_len = min(length(cancer), length(control))
    for i = 1:min_len
        println(out, cancer[i], "\t", control[i])
    end
    if length(cancer) > length(control)
        for i = (min_len+1):length(cancer)
            println(out, cancer[i])
        end
    else
        for i = (min_len+1):length(control)
            println(out, "\t", cancer[i])
        end
    end
    close(out)
end


println("Enter accession number:")
read_experiment(chomp(readline()), "data/dna-array-subset.dat")

#!/usr/bin/env julia
# Useful Functions and Methods
# https://teaching.healthtech.dtu.dk/22110/index.php/Useful_Functions_and_Methods

# Exercise 1
read_data() = try
    if length(ARGS) == 1
        filename = ARGS[1]
    else
        filename = "data/appendix1.txt"
    end
    open(_read_data, filename)
catch e
    msg = sprint(showerror, e)
    println(msg)
    exit(1)
end


function _read_data(f)
    SEQ_flag = false
    seq = ""
    muts = Array{Array{String, 1}, 1}()

    for line in eachline(f)
        FT_reg = match(r"^FT\s+(VARIANT|MUTAGEN)\s+([0-9]+)\s+\S+\s+(\S)\s*->\s*(\S)", line)

        if FT_reg !== nothing
            push!(muts, [FT_reg[x] for x in 1:4])
        end
        if startswith(line, "//")
            SEQ_flag = false
        end
        if SEQ_flag
            seq *= replace(strip(line), " " => "")
        end
        if startswith(line, "SQ")
            SEQ_flag = true
        end
    end

    println(muts)
    println(seq)

    out = open("mutations.fsa", "w")

    for mut in muts
        res_string = mut[2]
        res = parse(Int64, res_string)
        if mut[3] == string(seq[res])
            seq_mut = seq[begin:(res-1)] * mut[4] * seq[(res+1):end]
            println(out, ">", mut[1], mut[2], " ", mut[3], "->", mut[4])
            for i in 1:60:length(seq_mut)
                if i + 59 > length(seq_mut)
                    println(out, seq_mut[i:end])
                else
                    println(out, seq_mut[i:(i + 59)])
                end
            end
        else
            println(mut[3], " and ", seq[res], " did not match.")
        end
    end

    close(out)
end

read_data()


# Exercise 2
using DelimitedFiles

function readpoints(f)
    datapoints = Matrix{Float64}(undef, 6, 0)
    ids = Vector{String}()

    for ln in eachline(f)
        ln_split = split(ln, "\t")
        push!(ids, ln_split[1])
        datapoints = hcat(datapoints, [parse(Float64, x) for x in ln_split[2:end]])
    end

    ids, datapoints
end

function average(a)
    sum(a)/length(a)
end

function highest_lowest()
    translate_arr = readdlm("data/appendix5.txt", String)

    translate_dict = Dict(translate_arr[i, 1] => translate_arr[i, 3]
        for i in 1:size(translate_arr,1))

    exclude_arr = readlines("data/appendix4.txt")

    ids, data = readpoints("data/appendix3.txt")
    translated_exclude = Set([get(translate_dict, id, "") for id in exclude_arr])

    println("Excluded (translated): ", translated_exclude)

    avgs = Vector{Float64}(undef, length(ids))
    avgs_included = Int[]

    for i in eachindex(avgs)
        if ids[i] ∉ translated_exclude
            avgs[i] = average(data[:,i])
            push!(avgs_included, i)
        end
    end

    filtered_avgs = avgs[avgs_included]
    filtered_ids = ids[avgs_included]

    lowest_i = argmin(filtered_avgs)
    highest_i = argmax(filtered_avgs)

    println("Highest average: ", filtered_ids[highest_i], " : ", filtered_avgs[highest_i])
    println("Lowest average: ", filtered_ids[lowest_i], " : ", filtered_avgs[lowest_i])
end

highest_lowest()

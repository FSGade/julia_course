#!/usr/bin/env julia
# Lists/Sequences

# Exercise 1
let
    done = false
    words = String[]

    while !done
        println("Gimme some words. I will stop on 'STOP': ")
        word = chomp(readline())
        if word == "STOP"
            done = true
        else
            push!(words, word)
        end
    end
    open("words.txt", "w") do f
        for word in words
            println(f, word)
        end
    end
end


# Exercise 2

try
    words = readlines("words.txt")
    open("words.txt", "w") do f
        for word in reverse(sort(words))
            println(f, word)
        end
    end
catch e
    msg = sprint(showerror, e)
    println(msg)
end

# Exercise 3

try
    accessions = readlines("data/ex5.acc")
    sort!(accessions)
    open("clean.acc", "w") do f
        println(f, accessions[1])
        for i in 2:length(accessions)
            if accessions[i] != accessions[i-1]
                println(f, accessions[i])
            end
        end
    end
catch e
    msg = sprint(showerror, e)
    println(msg)
end


# Exercise 4
try
    accessions = readlines("data/ex5.acc")
    sort!(accessions)
    open("clean.acc", "w") do f
        i = 2
        while i <= length(accessions)
            if accessions[i] == accessions[i-1]
                popat!(accessions, i)
            else
                i += 1
            end
        end
        for acc in accessions
            println(f, acc)
        end
    end
catch e
    msg = sprint(showerror, e)
    println(msg)
end


# Exercise 5
function ask_for_acc()
    println("Search in clean.acc: ")
    chomp(readline())
end

function linear_search_clean()
    accessions = readlines("clean.acc")
    search = ask_for_acc()
    while search != "STOP"
        if search in accessions
            println(search*" was found in the file.")
        else
            println(search*" was NOT found in the file.")
        end
        search = ask_for_acc()
    end
    println("Done.")
end

linear_search_clean()

# Exercise 6
function binary_search_clean()
    accessions = readlines("clean.acc")
    search = ask_for_acc()
    while search != "STOP"
        minimum = 1
        maximum = length(accessions)
        inlist = false

        while minimum != maximum && !inlist
            mid = convert(Int64, floor((minimum + maximum) / 2))
            if accessions[mid] < search
                minimum = mid + 1
            elseif accessions[mid] > search
                maximum = mid
            else
                inlist = true
            end
        end
        inlist ? println(search*" was found in the file.") : println(search*" was NOT found in the file.")
        println()
        search = ask_for_acc()
    end
    println("Done.")
end

binary_search_clean()

# Exercise 7.a
function print_avg(nums)
    nums = [parse(Int64, num) for num in nums]
    println(sum(nums)/length(nums))
end

if length(ARGS) == 2
    print_avg(ARGS)
else
    nums = readlines()
    if length(nums) > 0
        print_avg(nums)
    else
        println("Input was not given correctly.")
    end
end



# Exercise 7.b

function count_negative(f)
    negative_count = 0
    for ln = eachline(f)
        negative_count += count("-", ln)
    end
    println(negative_count)
    negative_count
end

if length(ARGS) == 0
    println("Enter filename:")
    filename = chomp(readline())
    open(count_negative, filename)
elseif length(ARGS) == 1
    filename = ARGS[1]
    open(count_negative, filename)
else
    println("You need to supply only one filename.")
end


# Exercise 8
using DelimitedFiles

function print_help()
    print("Arguments given incorrectly")
    print("The program must be used as such:")
    print("julia $PROGRAM_FILE -fx,y file")
    exit(1)
end

function cut()
    if length(ARGS) == 2
        filename = ARGS[2]
        if ARGS[1][begin:2] == "-f"
            cols = [parse(Int64, x) for x in split(ARGS[1][3:end], ",")]
        else
            print_help()
        end
        data = readdlm(filename, '\t', Any, '\n')
        open("$filename.cut", "w") do io
           writedlm(io, data[:, cols])
       end
    else
        print_help()
    end
end

cut()

# Exercise 9-10

using DelimitedFiles

function print_help()
    print("Arguments given incorrectly")
    print("The program must be used as such:")
    print("julia $PROGRAM_FILE file")
    exit(1)
end

function sum_cols()
    if length(ARGS) == 1
        try
            data = readdlm(ARGS[1], '\t', Float64, '\n')
            num_points, num_cols = size(data)
            avg = Vector{Float64}(undef, num_cols)
            for i in eachindex(avg)
                avg[i] = sum(view(data, :, i))/num_points
            end
            println(avg)
        catch e
            msg = sprint(showerror, e)
            println(msg)
            exit(1)
        end
    else
        print_help()
    end
end

sum_cols()

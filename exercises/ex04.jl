#!/usr/bin/env julia
# Exceptions and Bug Handling
# Exercise 1-5

struct InvalidLengthError <: Exception
    val1::Int64
    val2::Int64
end
Base.showerror(io::IO, e::InvalidLengthError) = print(io, "The parsed length (", e.val1, ") does not equal the calculated length (",e.val2,")")


function parse_prot(f)
    seq = ""
    sp_id = ""
    sp_ac = ""
    aa_len = 0
    FLAG = false

    for line in eachline(f)
        key, value = line[begin:2], rstrip(line[6:end])
        if key == "ID"
            i = 1
            while value[i] != ' '
                i += 1
            end
            sp_id = string(value[begin:i])
            println("The ID is: ", sp_id)
        elseif key == "AC"
            sp_ac = string(value[begin:end-1])
            println("The accession numbers are: ", value[begin:end-1])
        end

        if key == "//"
            FLAG = false
        end
        if FLAG
            for char in value
                if char != ' '
                    seq *= char
                end
            end
        end
        if key == "SQ"
            FLAG = true
            aa_end = aa_start = 12
            while value[aa_end] != ' '
                aa_end += 1
            end
            aa_len = parse(Int64, value[aa_start:aa_end])
        end
    end
    println("Sequence: ",seq)

    if aa_len != length(seq)
        throw(InvalidLengthError(aa_len,length(seq)))
    else
        println("AA length: ", aa_len)
    end
    seq, aa_len, sp_ac, sp_id
end

read_prot(filename) = try
    open(parse_prot, filename)
catch e
    msg = sprint(showerror, e)
    println(msg)
end

println("Enter file containing SwissProt database entry to analyze: ")

seq, aa_len, sp_ac, sp_id = read_prot(chomp(readline()))

#Ex 6
function to_outfile(seq, aa_len, sp_ac, sp_id)
    outfile = open("sprot.fsa", "w")
    println(outfile,">sp|", sp_ac, "|", sp_id)

    for i in 1:60:length(seq)
        if i+59 > length(seq)
            println(outfile, seq[i:end])
        else
            println(outfile, seq[i:(i+59)])
        end
    end

    close(outfile)
end

to_outfile(seq, aa_len, sp_ac, sp_id)


# Exercise 7
function parse_dna(f)
    fasta_info = readline(f)
    seq = ""
    for line in eachline(f)
        seq *= chomp(line)
    end
    for i in 1:(length(seq)-2)
        println(seq[i:i+2])
        if seq[i:i+2] == "ATG"
            return seq, i
        end
    end
    seq, nothing
end

read_dna(filename = "data/dna.fsa") = try
    seq = open(parse_dna, filename)
catch e
    msg = sprint(showerror, e)
    println(msg)
end

println(read_dna())

# Exercise 8
function parse_dna(f)
    fasta_info = readline(f)
    seq = ""
    for line in eachline(f)
        seq *= chomp(line)
    end
    for i in 1:(length(seq)-2)
        println(seq[i:i+2])
        if seq[i:i+2] == "ATG"
            for j in i:3:(length(seq)-2)
                if seq[j:j+2] in ("TAA", "TAG", "TGA")
                    return seq, i, j
                end
            end
        end
    end
    seq, nothing, nothing
end

read_dna(filename = "data/dna.fsa") = try
    seq = open(parse_dna, filename)
catch e
    msg = sprint(showerror, e)
    println(msg)
end

println(read_dna())


# Exercise 9
function parse_orphans(match, f)
    match_count = 0
    for line in eachline(f)
        if line[1] != '>'
            i = 1
            while line[i] != '_'
                i += 1
            end
            if line[i+1:i+length(match)] == match
                match_count += 1
            end
        end
    end
    match_count
end

search_orphans(match, filename = "data/orphans.sp") = try
    seq = open(f->parse_orphans(match, f), filename)
catch e
    msg = sprint(showerror, e)
    println(msg)
end

println(search_orphans("MOUSE"))


# Exercise 10
function guess_number(minimum, maximum)
    ori_min = minimum
    ori_max = maximum

    gameFinished = false

    guess = convert(Int64, floor((minimum + maximum)/2))

    while !gameFinished
        if minimum == maximum # only 1 possible number left
            println("Your number is ", minimum, "!")
            gameFinished = true
        elseif minimum > maximum # the limits are weird, since you lied!
            println("You lied to me, you son of a gun! Bund and start over :)")
            minimum = ori_min
            maximum = ori_max
            #gameFinished = True
        else
            println("Is your number ", guess, "? If not, higher or lower?")
            choice = chomp(readline())
            if choice == "higher"
                minimum = guess + 1
            elseif choice == "lower"
                maximum = guess - 1
            elseif choice == "yes"
                println("Gotcha!")
                gameFinished = true
            end
        end
        guess = convert(Int64, floor((minimum + maximum)/2))
    end
end

guess_number(1, 10)

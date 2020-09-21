# Ex 3.1
println("First number: ")
num1 = parse(Int, chomp(readline()))
println("Second number: ")
num2 = parse(Int, chomp(readline()))

println(num1, " ", num2, " Average: ", Int(floor( (num1 + num2) / 2 )))

# Ex 3.2
# echo "54" > nums_file
# echo "72" >> nums_file
# cat nums_file | julia ex3.jl

# Ex 3.3

function count_negative(ln)
    count = 0
    match = findfirst(isequal('-'), ln)
    while match !== nothing
        count += 1
        match = findnext(isequal('-'), ln, match + 1)
    end
    count
end

num_neg = open("ex1.dat") do f
    num_neg = 0
    for ln in eachline(f)
        num_neg += count_negative(ln)
    end
    num_neg
end

# Ex 3.4
temperature = chomp(readline())
temp_num = parse(Int, temperature[1:end-1])
temp_scale = temperature[end]

if temp_scale == 'F'
    println((temp_num-32)/1.8, "C")
elseif temp_scale == 'C'
    println(temp_num*1.8 + 32, "F")
end

# Ex 3.5
function extract_acc(ln)
    last_index = findfirst(isequal(' '), ln)
    ln[2:last_index-1]
end

accs = open("/home/frederik/julia_course/data/orphans.sp") do f
    accs = String[]
    for ln in eachline(f)
        if ln[1] == '>'
            push!(accs, extract_acc(ln))
        end
    end
    accs
end

println(accs)

# Ex 3.6


# Ex 3.7-8
function read_dna(f)
    dna = ""
    for ln in eachline(f)
        dna = string(dna, rstrip(ln))
    end
    dna
end

in_dna = open(read_dna, "/home/frederik/julia_course/data/dna.dat")

function complement(str)
    trans = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C')
    trans_str = map(x -> trans[x], str)
end

function reverse_complement(str)
    reverse(complement(str))
end

println(in_dna)
println(reverse_complement(in_dna))

# Ex 3.9

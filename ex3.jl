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
function to_outfile(seq)
    outfile = open("revdna.dat", "w")

    for i in 1:60:length(seq)
        if i+59 > length(seq)
            println(outfile, seq[i:end])
        else
            println(outfile, seq[i:(i+59)])
        end
    end

    close(outfile)
end

to_outfile(reverse_complement(in_dna))

# Ex 3.10
function read_fasta(f)
    headlines = String[]
    seqs = String[]
    seq = ""
    for line = eachline(f)
        ln = chomp(line)
        if ln[1] == '>'
            if seq != ""
                push!(seqs, seq)
                seq = ""
            end
            push!(headlines, ln)
        else
            seq *= ln
        end
    end
    push!(seqs, seq)
    headlines, seqs
end

open("data/revdna.fsa", "w") do outfile
    headlines, seqs = open(read_fasta, "data/dna.fsa")
    for (headline, seq) in zip(headlines, seqs)
        println(outfile, "$headline ComplementStrand")
        revseq = reverse_complement(seq)
        for i in 1:60:length(revseq)
            if i + 59 > length(revseq)
                println(outfile, revseq[i:end])
            else
                println(outfile, revseq[i:(i + 59)])
            end
        end
    end
end

# Exercise 11
function count_ATCG(seq)
    AT = 0
    GC = 0
    for char in seq
        if char in ('A', 'T')
            AT += 1
        elseif char in ('G', 'C')
            GC += 1
        end
    end
    return string("AT/GC content: ", AT/GC)
end


headlines, seqs = open(read_fasta, "data/dna.fsa")
println(count_ATCG(seqs[1]))

# Exercise 12
function bullseye(dimensions = (40, 40), c1 = "*", c2 = "-", c3 = " ")
    xdim, ydim = dimensions
    for x in 1:xdim
        for y in 1:ydim
            _r = sqrt((x-xdim/2)^2+(y-ydim/2)^2)
            if _r < 5
                print(c1)
            elseif _r < 20
                print(c2)
            else
                print(c3)
            end
        end
        println()
    end
end

bullseye()

#!/usr/bin/env julia
# Pattern Matching and Regular Expressions
# https://teaching.healthtech.dtu.dk/22110/index.php/Pattern_Matching_and_Regular_Expressions

# Exercise 1
println("Check if is number:")
num = chomp(readline())
if match(r"^-?\d+(\.\d+)?$", num) !== nothing
    println("Yes, it is a number")
else
    println("No, it is not a number")
end

# Exercise 2
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
            id_reg = match(r"^(\w+)\s+", value)
            sp_id = string(id_reg[1])
            println("The ID is: ", sp_id)
        elseif key == "AC"
            ac_reg = match(r"^(\w+);", value)
            sp_ac = string(ac_reg[1])
            println("The accession number is: ", sp_ac)
        end

        if key == "//"
            FLAG = false
        end
        if FLAG
            seq_parts_reg = findall(r"(\w+)", value)
            seq_parts = [value[x] for x in seq_parts_reg]
            seq *= join(seq_parts)
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

# Exercise 3
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
    try
        headlines, seqs = open(read_fasta, "data/dna.fsa")
        for (headline, seq) in zip(headlines, seqs)
            if length(findfirst(r"[ATCG]+", seq)) != length(seq)
                throw(ArgumentError("The sequence does not entirely contain ATCG."))
            end
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
    catch err
        msg = sprint(showerror, err)
        print(msg)
        exit(1)
    end
end

# Exercise 4-9
read_data() = try
    if length(ARGS) == 1
        filename = ARGS[1]
    else
        filename = "data/data1.gb"
    end
    open(_read_data, filename)
catch e
    msg = sprint(showerror, e)
    println(msg)
    exit(1)
end

function _read_data(f)
    MED = String[]
    SEQ_flag = AA_flag = CDS_flag = false
    SEQ = AA = CDS = ""
    CDS_range = ""

    for line in eachline(f)
        if startswith(line, "//")
            SEQ_flag = false
        end
        if SEQ_flag
            SEQ_re = findall(r"([atcg]+)", line)
            SEQ_parts = [line[x] for x in SEQ_re]
            SEQ *= join(SEQ_parts)
        else
            ACC_re = match(r"^ACCESSION\s+(\S+)$", line)
            DEF_re = match(r"^DEFINITION\s+(.+)$", line)
            ORG_re = match(r"^\s+ORGANISM\s+(.+)$", line)
            MED_re = match(r"^\s+MEDLINE\s+(\S+)$", line)

            if ACC_re !== nothing
                println("Accession number: ", ACC_re[1])
            elseif DEF_re !== nothing
                println("Definition: ", DEF_re[1])
            elseif ORG_re !== nothing
                println("Organism: ", ORG_re[1])
            elseif MED_re !== nothing
                push!(MED, MED_re[1])
            end
        end

        if startswith(line, "ORIGIN")
            SEQ_flag = true
        end

        CDS_re = match(r"^\s+CDS\s+(\S+)$", line)
        AA_re = match(r"^\s+/translation=(\S+)$", line)

        if CDS_flag
            CDS_data = strip(line)
            if startswith(CDS_data, "/")
                CDS_flag = false
            else
                CDS_range *= CDS_data
            end
        end
        if CDS_re !== nothing
            CDS_range *= CDS_re[1]
            CDS_flag = true
        end

        if AA_flag
            AA_data = match(r"^\s+(\S+)$", line)
            if AA_data !== nothing
                AA *= AA_data[1]
            else
                AA_flag = false
            end
        end
        if AA_re !== nothing
            AA *= AA_re[1]
            AA_flag = true
        end
    end

    CDS_list = [split(CDS_range[x], "..") for x in findall(r"(\d+\.\.\d+)", CDS_range)]

    for x in CDS_list
        CDS *= SEQ[parse(Int64, x[1]):parse(Int64, x[2])]
    end

    println("MEDLINE article number(s): ", join(MED, ", "))
    println()
    println("Sequence: ", SEQ)
    println()
    println("Coding sequence (CDS): ", CDS)
    println()
    println("Amino acid sequence: ", AA)
end

read_data()

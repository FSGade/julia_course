#!/usr/bin/env julia
# Sets and Dictionaries
# https://teaching.healthtech.dtu.dk/22110/index.php/Sets_and_Dictionaries

const codons = Dict("TTT"=>"F", "TTC"=>"F", "TTA"=>"L", "TTG"=>"L",
    "TCT"=>"S", "TCC"=>"S", "TCA"=>"S", "TCG"=>"S",
    "TAT"=>"Y", "TAC"=>"Y", "TAA"=>"", "TAG"=>"",
    "TGT"=>"C", "TGC"=>"C", "TGA"=>"", "TGG"=>"W",
    "CTT"=>"L", "CTC"=>"L", "CTA"=>"L", "CTG"=>"L",
    "CCT"=>"P", "CCC"=>"P", "CCA"=>"P", "CCG"=>"P",
    "CAT"=>"H", "CAC"=>"H", "CAA"=>"Q", "CAG"=>"Q",
    "CGT"=>"R", "CGC"=>"R", "CGA"=>"R", "CGG"=>"R",
    "ATT"=>"I", "ATC"=>"I", "ATA"=>"I", "ATG"=>"M",
    "ACT"=>"T", "ACC"=>"T", "ACA"=>"T", "ACG"=>"T",
    "AAT"=>"N", "AAC"=>"N", "AAA"=>"K", "AAG"=>"K",
    "AGT"=>"S", "AGC"=>"S", "AGA"=>"R", "AGG"=>"R",
    "GTT"=>"V", "GTC"=>"V", "GTA"=>"V", "GTG"=>"V",
    "GCT"=>"A", "GCC"=>"A", "GCA"=>"A", "GCG"=>"A",
    "GAT"=>"D", "GAC"=>"D", "GAA"=>"E", "GAG"=>"E",
    "GGT"=>"G", "GGC"=>"G", "GGA"=>"G", "GGG"=>"G")

function translate_dna(dna)
    if length(dna) % 3 != 0
        println("Whoops")
    end
    peptide = ""
    for i = 1:3:length(dna)
        peptide *= codons[dna[i:i+2]]
    end
    peptide
end

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
            seq *= translate_dna(ln)
        end
    end
    push!(seqs, seq)
    headlines, seqs
end

open("data/aa7.fsa", "w") do outfile
    headlines, seqs = open(read_fasta, "data/dna7.fsa")
    for (headline, seq) in zip(headlines, seqs)
        println(outfile, "$headline Amino Acid Sequence")
        for i in 1:60:length(seq)
            if i + 59 > length(seq)
                println(outfile, seq[i:end])
            else
                println(outfile, seq[i:(i + 59)])
            end
        end
    end
end

# Ex 3
function read_as_set(f, colnum)
    extract_line = if colnum === 0
        ef1(ln) = string(chomp(ln))
    else
        ef2(ln) = string(split(ln, "\t")[colnum])
    end
    outset = Set{String}()
    for line = eachline(f)
        push!(outset, extract_line(line))
    end
    outset
end

function read_set_diff(f1, f2)
    set1 = read_as_set(f1, 0)
    set2 = read_as_set(f2, 2)
    setdiff(set1,set2)
end

println(read_set_diff("data/start10.dat", "data/res10.dat"))

# Ex 4
function read_as_count_dict(f)
    outdict = Dict{String, Int64}()
    for line = eachline(f)
        x = chomp(line)
        outdict[x] = get(outdict, x, 0) + 1
    end
    outdict
end

open("data/ex5.acc.count", "w") do f
    for (id, count) in sort(collect(open(read_as_count_dict, "data/ex5.acc")),
            by=x->x[2], rev=true)
        println(f, id, " ", count)
    end
end

# Ex 5
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

    CDS_dict = Dict{String, Int64}()
    for i = 1:3:length(CDS)
        codon = CDS[i:(i+2)]
        CDS_dict[codon] = get(CDS_dict, codon, 0) + 1
    end


    println("MEDLINE article number(s): ", join(MED, ", "))
    println()
    println("Sequence: ", SEQ)
    println()
    println("Coding sequence (CDS): ", CDS)
    println()
    println("Amino acid sequence: ", AA)
    println()
    println("Codons: ")
    for (codon, count) in sort(collect(CDS_dict), by=x->x[2], rev=true)
        println(codon, ": ", count)
    end
end

read_data()


# Ex 6

read_authors() = try
    if length(ARGS) == 1
        filename = ARGS[1]
    else
        filename = "data/data1.gb"
    end
    open(_read_authors, "data/data1.gb")
catch e
    msg = sprint(showerror, e)
    println(msg)
    #exit(1)
end

function _read_authors(f)
    MED = String[]
    MED_flag = false

    for line in eachline(f)
        MED_reg = match(r"^\s+AUTHORS\s+(.*)$", line)

        if MED_flag
            MED_extra_reg = match(r"^\s+[A-Z]+\s", line)
            if MED_extra_reg === nothing
                MED[end] *= strip(line)
            else
                MED_flag = false
            end
        end

        if MED_reg !== nothing
            push!(MED, MED_reg[1])
            MED_flag = true
        end
    end

    MED = [replace(x, " and" => ", ") for x in MED]

    authors_arr = String[]
    for names in MED
        append!(authors_arr, [strip(name) for name in split(names, ", ")])
    end

    authors = Set(authors_arr)

    println(MED)
    println(authors)
    authors
end

read_authors()

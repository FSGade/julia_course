#=
# Ex 2.1
open("ex1.acc") do file
    for ln in eachline(file)
        println(ln)
    end
end

# Ex 2.2
println("Enter a file name:")
fn = chomp(readline())

open(fn) do f
    for ln in eachline(f)
        println(ln)
    end
end

# Ex 2.3
println("Enter a file name:")
fn = chomp(readline())

nl = open(fn) do f
    i = 0
    for ln in eachline(f)
        i += 1
    end
    i
end

println(nl)
=#

#Ex 2.4-9
for i = 1:3
    avg, sum, lc, neg, pos, zer, max_num, min_num = open("ex1.dat.$i") do f
        lc, sum, neg, pos, zer = 0, 0, 0, 0, 0
        max_num, min_num = -Inf, Inf
        for ln in eachline(f)
            ln_num = parse(Float64, chomp(ln))
            (ln_num > 0) ? pos += 1 : (ln_num < 0) ? neg += 1 : zer += 1
            (ln_num > max_num) && (max_num = ln_num)
            (ln_num < min_num) && (min_num = ln_num)
            lc += 1
            sum += ln_num
        end
        sum / lc, sum, lc, neg, pos, zer, max_num, min_num
    end
    println("Column ", i, " - Average: ", avg, ", Sum: ", sum,
    ", Line count: ", lc)
    println("Negative: ", neg, ", Positive: ", pos, ", Zeros: ", zer)
    println("Maximum number: ", max_num, ", Minimum number: ", min_num)
    println()
end

#=
avg = open("1to9") do f
    i, sum = 0, 0
    for ln in eachline(f)
        i += 1
        sum += parse(Float64, chomp(ln))
    end
    sum / i
end
println(avg)
=#

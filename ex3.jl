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
temp_scale = string(temperature[end])

if temp_scale == "F"
    println((temp_num-32)/1.8, "C")
elseif temp_scale == "C"
    println(temp_num*1.8 + 32, "F")
end

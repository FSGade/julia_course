
# Ex 1.1
println("Hello World")

# Ex 1.2
for i in 1:10
    println("Hello World")
end

# Ex 1.3
for i in 1:10
    println(i)
end

# Ex 1.4
println("What is your name?")
x = chomp(readline())
if x == "Frederik"
    println("That is a good name, ", x, ". Welcome!")
else
    println("Welcome ", x, "!")
end

# Ex 1.5
println("Enter a number")
x = parse(Float64, chomp(readline()))
println("Enter another number")
y = parse(Float64, chomp(readline()))
println(x, " + ", y, " = ", x + y)

# Ex 1.6
println("Enter a number")
x = parse(Float64, chomp(readline()))
println("Enter another number")
y = parse(Float64, chomp(readline()))

println("Enter an operator (+, -, *, /)")
a = chomp(readline())
f = getfield(Main, Symbol(a))
println(x, " ", a, " ", y, " = ", f(x,y))

# Ex 1.7
println("Enter a number")
x = parse(Int64, chomp(readline()))
println("Enter another number")
y = parse(Int64, chomp(readline()))

for i in x:y
    println(i)
end

# Ex 1.8
println("Enter a number")
x = parse(Int64, chomp(readline()))
println("Enter another number")
y = parse(Int64, chomp(readline()))

if x > y
    x, y = y, x
end

for i in x:y
    println(i)
end

# Ex 1.9
function gt_previous(val, list)
    for elem in list
        if val < elem
            return false
        end
    end
    true
end

println("Enter a number")
x = parse(Int64, chomp(readline()))
num_list = []

while gt_previous(x, num_list)
    push!(num_list, x)
    println("Enter a number")
    global x = parse(Int64, chomp(readline()))
end


# Ex 1.10
println("Enter a non-negative number")
x = parse(Int64, chomp(readline()))

function fact(n::Int64)::Int64
    if n < 0
        error("Negative values are not allowed for factorial.")
    end
    res = 1
    for i in 1:n
        res *= i
    end
    return res
end
println(fact(x))


println("Enter a non-negative number")
x = parse(Int64, chomp(readline()))

fact(n) = (n > 0) ? *(1:n...) :
            (n < 0) ? error("Only non-negative vals allowed.") : 1

println(fact(x))


# Ex 1.11
println("Enter a number")
x = parse(Int64, chomp(readline()))

#if x > 0
#    r = 1:x
#else
#    r = x:-1
#end

#sum = +(r...)

sum = +(1:abs(x)...) * sign(x)

println(sum)

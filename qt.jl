#!/usr/bin/env julia
# PRELIMINARY QT CLUSTERING IMPLEMENTATION
# NO OPTIMISATIONS (ONLY "DIAMETER CACHE")
# AUTHOR: FREDERIK GADE

@show ARGS
println(ARGS)

function dp(lnsplit, offset)
    [parse(Float64, x) for x in lnsplit[offset:end]]
end

#TODO: Name list
function read_points(f)
    i = 1
    firstline = split(readline(f), "\t")
    if tryparse(Float64, firstline[1]) === nothing
        names = [firstline[1]]
        dim = length(firstline) - 1
        offset = 2
        datapoints = dp(firstline, offset)
    else
        names = [string("Point", i)]
        dim = length(firstline)
        offset = 1
        datapoints = dp(firstline, offset)
    end
    for ln in eachline(f)
        datapoints = hcat(datapoints, dp(split(ln, "\t"), offset))
    end
    names, datapoints
end

function squared!(a)
    @simd for i in eachindex(a)
        @inbounds a[i] = a[i]^2
    end
    a
end

function dist(x, y, d)
    r = Vector{Float64}(undef, d)
    @simd for i = 1:d
        @inbounds r[i] = y[i] - x[i]
    end
    squared!(r) |> (sqrt ∘ sum)
end

function get_dists(datapoints)
    dim, num_points = size(datapoints)
    diameter = -Inf
    dists = zeros(eltype(datapoints), num_points, num_points)
    for i = 1:num_points
        for j = i+1:num_points
            d = dist(datapoints[:,i], datapoints[:,j], dim)
            @inbounds dists[i,j] = dists[j,i] = d
            if d > diameter
                diameter = d
            end
        end
    end
    dists, diameter
end

function get_threshold(percentage, diameter)
    diameter*percentage/100
end

function get_neighbours(dists, threshold, num_points)
    neighbours = [Int[] for i = 1:num_points]
    for i = 1:num_points
        @simd for j = i+1:num_points
            if dists[i, j] <= threshold
                push!(neighbours[i], j)
                push!(neighbours[j], i)
            end
        end
    end
    neighbours
end

function generate_candidate(seed, dists, neighbours, threshold)
    candidate_cluster = [seed]
    last_added = seed
    seed_neighbours = copy(neighbours[seed])
    diameter_cache = zeros(length(seed_neighbours))
    while length(seed_neighbours) > 0
        for index in 1:length(seed_neighbours)
            d = dists[last_added, seed_neighbours[index]]
            if diameter_cache[index] < d
                diameter_cache[index] = d
            end
        end

        candidate_point = argmin(diameter_cache)
        candidate_point_val = seed_neighbours[candidate_point]

        if diameter_cache[candidate_point] < threshold
            push!(candidate_cluster, candidate_point_val)
            last_added = candidate_point_val
            popat!(seed_neighbours, candidate_point)
            popat!(diameter_cache, candidate_point)
        else
            return candidate_cluster
        end
    end
    return candidate_cluster
end

function remove_clustered(remove, arr)
    filter!(x-> x ∉ remove, arr)
end

function update_data(best, unclustered, neighbours)
    remove_clustered(best, unclustered)
    remove_clustered(best, neighbours)
    for i in eachindex(neighbours)
        remove_clustered(best, neighbours[i])
    end
    unclustered, neighbours
end

function update_output(output, best_cluster)
    push!(output, sort(map(x -> x-1, best_cluster)))
end

function QT(filename, percentage)
    names, datapoints = open(read_points, filename)
    num_points = size(datapoints, 2)
    dists, diameter = get_dists(datapoints)
    threshold = get_threshold(percentage, diameter)
    neighbours = get_neighbours(dists, threshold, num_points)

    unclustered = collect(1:num_points)
    output = []

    while length(unclustered) != 0
        best_cluster_size = 0
        best_cluster = []

        for seed in unclustered
            candidate = generate_candidate(seed, dists, neighbours, threshold)
            candidate_size = length(candidate)
            if candidate_size > best_cluster_size
                best_cluster = candidate
                best_cluster_size = candidate_size
            end
        end

        update_output(output, best_cluster)
        update_data(best_cluster, unclustered, neighbours)
    end
    output
end

println("---------------------------")
@time QT("/home/frederik/Desktop/julia_course/data/point1000.lst", 30) |> println

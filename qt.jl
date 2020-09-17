#!/usr/bin/env julia
# PRELIMINARY QT CLUSTERING IMPLEMENTATION
# SAVE CANDIDATE CLUSTERS & "DIAMETER CACHE" OPTIMISATIONS
# AUTHOR: FREDERIK GADE

#@show ARGS

###
### HELPER FUNCTIONS
###
function get_args()
    arg_len = length(ARGS)
    if arg_len == 2
        filename = ARGS[1]
        threshold = ARGS[2]
    elseif arg_len == 1
        if ARGS[1] in ("-h","--help")
            println("usage: julia "*PROGRAM_FILE* " inputfile threshold")
            println()
            println("""QT (Quality Threshold) Clustering is an algorithm that groups multi-dimensional
            vectors into high quality clusters. Quality is ensured by finding large cluster
            whose diameter does not exceed a given user-defined diameter threshold.
            This method prevents dissimilar vectors from being forced under the same
            cluster and ensures that only good quality clusters will be formed.""")
            println()
            println("arguments:")
            println("""    inputfile\tInput file name (tab-separated list of floats w/ optional name
            \t\tin first column)""")
            println("    threshold\tQuality threshold of cluster diameter")
            println()
        else
            println("Only one argument was supplied. Use flag -h/--help for help.")
        end
        exit(1)
    elseif arg_len == 0
        println("No arguments supplied. Use flag -h/--help for help.")
        exit(1)
    else
        println("Too many arguments supplied (only 2 allowed). Use flag -h/--help for help.")
        exit(1)
    end
    return filename, parse(Int64, threshold)
end

function get_threshold(percentage::Int64, diameter::Float64)::Float64
    diameter * percentage / 100
end

function tab_split(ln)
    split(ln, "\t")
end

function dp(lnsplit::Array{SubString{String},1},
        offset::Int64, i::Int64)
    datapoint = [parse(Float64, x) for x in lnsplit[offset:end]]
    if offset == 1
        return string("Point", lpad(i, 4, "0")), datapoint
    else
        return string(lnsplit[1]), datapoint
    end
end

function remove_clustered!(arr::Array{Int64,1}, remove::Array{Int64,1})
    filter!(x-> x ∉ remove, arr)
end
function remove_clustered!(arr::Array{Array{Int64,1},1},
        remove::Array{Int64,1})
    filter!(x-> x ∉ remove, arr)
end

function sq_diff!(r,x,y)
    @inbounds @simd for i = eachindex(r)
        @inbounds r[i] = (y[i] - x[i])^2
    end
end

#TODO: LISPify with exprs to simplify to dim might optimise?
function dist(x::Array{Float64,1}, y::Array{Float64,1}, d::Int64)::Float64
    #ntuple(i->(y[i] - x[i])^2, d) |> (sqrt ∘ sum)
    r = Vector{Float64}(undef, d)
    sq_diff!(r,x,y)
    sqrt(sum(r))
end

###
### PREPROCESSING FUNCTIONS
###

function validate_input(filename, percentage)
    println()
end

#TODO: Name list
function read_points(f)
    i = 1
    firstline = tab_split(readline(f))
    if tryparse(Float64, firstline[1]) === nothing
        dim = length(firstline) - 1
        offset = 2
    else
        dim = length(firstline)
        offset = 1
    end

    name, datapoint = dp(firstline, offset, i)
    names = [name]
    datapoints = Matrix{Float64}(undef, dim, 0)
    datapoints = hcat(datapoints, datapoint)

    for ln in eachline(f)
        i += 1
        name, datapoint = dp(tab_split(ln), offset, i)
        datapoints = hcat(datapoints, datapoint)
        push!(names, name)
    end

    names, datapoints
end

function get_dists(datapoints)
    dim, num_points = size(datapoints)
    diameter = -Inf
    dists = zeros(eltype(datapoints), num_points, num_points)
    @inbounds for i = 1:num_points
        @inbounds for j = i+1:num_points
            @inbounds d = dist(datapoints[:,i], datapoints[:,j], dim)
            @inbounds dists[i,j] = dists[j,i] = d
            if d > diameter
                diameter = d
            end
        end
    end
    dists, diameter
end

function get_neighbours(dists::Array{Float64,2}, threshold::Float64,
        num_points::Int64)
    neighbours = [Int[] for i = 1:num_points]
    for i = 1:num_points
        @simd for j = i+1:num_points
            @inbounds if dists[i, j] <= threshold
                push!(neighbours[i], j)
                push!(neighbours[j], i)
            end
        end
    end
    neighbours
end

###
### CORE ALGORITHM FUNCTIONS
###
function update_diameter_cache!(diameter_cache, dists, seed_neighbours,
        last_added, threshold)
    @inbounds d = dists[seed_neighbours, last_added]
    min_dist = threshold +1
    min_i = 1
    for (i, (x, y)) in enumerate(zip(diameter_cache, d))
        diameter_cache[i] = x > y ? x : y
        (min_i, min_dist) = diameter_cache[i] < min_dist ?
            (i, diameter_cache[i]) : (min_i, min_dist)
    end
    min_i
end

function update_data!(unclustered::Array{Int64,1},
        neighbours::Array{Array{Int64,1},1},
        candidate_clusters::Array{Array{Int64,1},1}, best::Array{Int64,1})
    remove_clustered!(unclustered, best)
    remove_clustered!(neighbours, best)
    remove_clustered!(candidate_clusters, best)
    for i in eachindex(neighbours)
        remove_clustered!(neighbours[i], best)
        if length(intersect(candidate_clusters[i], best)) != 0
            candidate_clusters[i] = Int[]
        end
    end
    nothing
end

function generate_candidate(seed, dists, neighbours::Array{Array{Int64,1},1},
        threshold::Float64)
    candidate_cluster = [seed]
    last_added = seed
    seed_neighbours = copy(neighbours[seed])
    diameter_cache = zeros(Float64, length(seed_neighbours))
    candidate_diameter = zero(Float64)

    while length(seed_neighbours) > 0
        candidate_point = update_diameter_cache!(diameter_cache, dists,
            seed_neighbours, last_added, threshold)
        candidate_point_val = seed_neighbours[candidate_point]

        if diameter_cache[candidate_point] < threshold
            candidate_diameter = diameter_cache[candidate_point]
            @inbounds push!(candidate_cluster, candidate_point_val)
            @inbounds popat!(seed_neighbours, candidate_point)
            @inbounds popat!(diameter_cache, candidate_point)
            last_added = candidate_point_val
        else
            return candidate_cluster, candidate_diameter
        end
    end

    if candidate_diameter == 0
        return candidate_cluster, Inf
    end

    return candidate_cluster, candidate_diameter
end

function write_output(outfile, i, names, datapoints)
    println(outfile, "-> Cluster ", i)
    for j in eachindex(names)
        print(outfile, names[j])
        for point in datapoints[:, j]
            print(outfile, "\t", point)
        end
        println(outfile)
    end
end

###
### MAIN FUNCTION
###
function QT(filename::String, percentage::Int64)
    validate_input(filename, percentage)
    names, datapoints = open(read_points, filename)
    num_points = size(datapoints, 2)
    dists, diameter = get_dists(datapoints)
    threshold = get_threshold(percentage, diameter)
    neighbours = get_neighbours(dists, threshold, num_points)

    unclustered = collect(1:num_points)
    candidate_clusters = [Int[] for i = 1:num_points]
    candidate_diameters = Vector{Float64}(undef, num_points)

    cluster_num = 1
    outfile = open(filename*".out", "w")

    while length(unclustered) != 0
        best_cluster_size = 0
        best_diameter = 0
        best_cluster = Int[]

        for seed in unclustered
            if length(candidate_clusters[seed]) != 0
                candidate, candidate_diameter =
                    candidate_clusters[seed], candidate_diameters[seed]
                candidate_size = length(candidate)
            else
                candidate, candidate_diameter =
                    generate_candidate(seed, dists, neighbours, threshold)
                candidate_clusters[seed] = candidate
                candidate_diameters[seed] = candidate_diameter
                candidate_size = length(candidate)
            end

            if candidate_size > best_cluster_size ||
                    (candidate_size == best_cluster_size
                    && candidate_diameter < best_diameter)
                best_cluster = candidate
                best_cluster_size = candidate_size
                best_diameter = candidate_diameter
            end
        end

        write_output(outfile, cluster_num,
            names[best_cluster], datapoints[:,best_cluster])
        cluster_num += 1

        update_data!(unclustered, neighbours, candidate_clusters, best_cluster)
    end

    close(outfile)
end

#@profiler @time QT(get_args()...)
@time QT(get_args()...)

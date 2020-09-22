#!/usr/bin/env julia
# PRELIMINARY QT CLUSTERING IMPLEMENTATION
# SAVE CANDIDATE CLUSTERS & "DIAMETER CACHE" OPTIMISATIONS
# AUTHOR: FREDERIK GADE

using LinearAlgebra: dot, mul!

###
### HELPER FUNCTIONS
###
function get_args()
    arg_len = length(ARGS)
    if arg_len == 2
        filename = ARGS[1]
        _threshold = ARGS[2]
    elseif arg_len == 1
        if ARGS[1] in ("-h", "--help")
            println("usage: julia ", PROGRAM_FILE, " inputfile threshold")
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
    return filename, _threshold
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

###
### PREPROCESSING FUNCTIONS
###

function validate_input(filename, _threshold)
    percentage = false
    lc = 0
    try
        lc = countlines(filename)
    catch e
        msg = sprint(showerror, e)
        println(msg)
    end
    threshold = tryparse(Float64, _threshold)
    if threshold === nothing
        if _threshold[end] == '%'
            percentage = true
            _threshold = chop(_threshold)
            threshold = tryparse(Float64, _threshold)
            (threshold === nothing) && throw(ArgumentError())
        end
    end
    threshold, percentage, lc
end

#TODO: Name list
function read_points(f, lc)
    i = 1
    firstline = tab_split(readline(f))
    if tryparse(Float64, firstline[1]) === nothing
        dim = length(firstline) - 1
        offset = 2
    else
        dim = length(firstline)
        offset = 1
    end

    datapoints = Matrix{Float64}(undef, dim, lc)
    names = Vector{String}(undef, lc)

    names[i], datapoints[:, i] = dp(firstline, offset, i)
    for ln in eachline(f)
        i += 1
        names[i], datapoints[:, i] = dp(tab_split(ln), offset, i)
    end

    names, datapoints
end

function sumsq_percol(a::Array{Float64,2})
    n = size(a, 2)
    r = Vector{Float64}(undef, n)
    @simd for j in 1:n
        aj = view(a, :, j)
        r[j] = dot(aj, aj)
    end
    return r
end

function _pairwise!(r::Array{Float64,2}, a::Array{Float64,2})
    m, n = size(a)
    mul!(r, a', a)
    sa2 = sumsq_percol(a)
    @inbounds for j = 1:n
        for i = 1:(j - 1)
            r[i, j] = r[j, i]
        end
        r[j, j] = 0
        sa2j = sa2[j]
        @simd for i = (j + 1):n
            r[i, j] = sqrt(max(sa2[i] + sa2j - 2 * r[i, j], 0))
        end
    end
    r
end

function get_dists(datapoints)
    n = size(datapoints, 2)
    r = Matrix{Float64}(undef, n, n)
    _pairwise!(r, datapoints)
    diameter = maximum(r)
    r, diameter
end

function get_neighbours(dists::Array{Float64,2}, threshold::Float64,
        num_points::Int64)
    neighbours = [Int[] for i = 1:num_points]
    for i = 1:num_points
        @simd for j = i+1:num_points
            @inbounds if dists[j, i] <= threshold
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
    filter!(x-> x ∉ best, unclustered)
    for i in unclustered
        neighbours[i] = setdiff(neighbours[i], best)
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

function get_best_cluster(unclustered, dists, neighbours, threshold,
        candidate_clusters, candidate_diameters)
    best_cluster_size = 0
    best_diameter = 0
    best_cluster_size = zero(Int64)
    best_diameter = zero(Float64)
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
    best_cluster
end

###
### MAIN FUNCTION
###
function QT(filename::String, _threshold::String)
    threshold, percentage, lc = validate_input(filename, _threshold)
    names, datapoints = open(f->read_points(f, lc), filename)
    num_points = size(datapoints, 2)
    dists, diameter = get_dists(datapoints)
    if percentage
        threshold = threshold/100 * diameter
    end
    neighbours = get_neighbours(dists, threshold, num_points)

    unclustered = collect(1:num_points)
    candidate_clusters = [Int[] for i = 1:num_points]
    candidate_diameters = Vector{Float64}(undef, num_points)

    cluster_num = 1
    outfile = open(filename*".out", "w")

    while length(unclustered) > 0
        best_cluster = get_best_cluster(unclustered, dists, neighbours,
            threshold, candidate_clusters, candidate_diameters)

        #sorted_cluster = sort(best_cluster)
        sorted_cluster = best_cluster

        write_output(outfile, cluster_num,
            names[sorted_cluster], datapoints[:,sorted_cluster])
        cluster_num += 1

        update_data!(unclustered, neighbours, candidate_clusters, best_cluster)
    end

    close(outfile)
end

#@profiler @time QT(get_args()...)
@time QT(get_args()...)

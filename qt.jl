#!/usr/bin/env julia
# QT CLUSTERING IMPLEMENTATION
# AUTHOR: FREDERIK GADE

using LinearAlgebra: dot, mul!

const SAVE_CLUSTER_VALS = (2, 3, 6, 9, 11, 23, 47, 106, 235, 551, 699)
const HELP_MESSAGE = """usage: julia $PROGRAM_FILE inputfile threshold

QT (Quality Threshold) Clustering is an algorithm that groups multi-dimensional
vectors into high quality clusters. Quality is ensured by finding large cluster
whose diameter does not exceed a given user-defined diameter threshold.
This method prevents dissimilar vectors from being forced under the same
cluster and ensures that only good quality clusters will be formed.

arguments:
    inputfile\tInput file name (tab-separated list of floats w/ optional name
\t\tin first column)
    threshold\tQuality threshold of cluster diameter
"""
const NO_ARGS = "No arguments supplied. Use -h/--help for help."
const ONE_ARG = "Only one argument was supplied. Use -h/--help for help."
const TOO_MANY_ARGS = "Too many arguments supplied. Use -h/--help for help."

##############################
### INPUT/OUTPUT FUNCTIONS ###
##############################
function tab_split(ln)
    split(ln, "\t")
end

function dp(lnsplit::Array{SubString{String},1},
        offset::Int64, i::Int64)
    datapoint = [parse(Float64, x) for x in lnsplit[offset:end]]
    if offset === 1
        return string("Point", lpad(i, 4, "0")), datapoint
    else
        return string(lnsplit[1]), datapoint
    end
end

function validate_input(filename, _threshold)
    percentage = false
    num_points = countlines(filename)
    threshold = tryparse(Float64, _threshold)
    if threshold === nothing
        if _threshold[end] == '%'
            percentage = true
            _threshold = chop(_threshold)
            threshold = tryparse(Float64, _threshold)
            (threshold === nothing) &&
                throw(ArgumentError("Threshold has to be a number/percentage"))
            if threshold < 0
                throw(ArgumentError("Threshold has to be positive"))
            elseif threshold > 100
                @warn "Percentage given is over 100%."
            end
        end
    end
    threshold, percentage, num_points
end

function read_points(f, num_points)
    i = 1
    firstline = tab_split(readline(f))
    if tryparse(Float64, firstline[1]) === nothing
        dim = length(firstline) - 1
        offset = 2
    else
        dim = length(firstline)
        offset = 1
    end

    datapoints = Matrix{Float64}(undef, dim, num_points)
    names = Vector{String}(undef, num_points)

    names[i], datapoints[:, i] = dp(firstline, offset, i)
    for ln in eachline(f)
        i += 1
        names[i], datapoints[:, i] = dp(tab_split(ln), offset, i)
    end

    names, datapoints
end

function write_output(outfile::IO, i::Int64, names::Array{String, 1},
    datapoints::Array{Float64,2})
    println(outfile, "-> Cluster ", i)
    for j in eachindex(names)
        println(outfile, join(vcat(names[j], datapoints[:, j]), "\t"))
    end
end

###############################
### PREPROCESSING FUNCTIONS ###
###############################
function sumsq_percol(a::Array{Float64,2})
    n = size(a, 2)
    r = Vector{Float64}(undef, n)
    @simd for j in 1:n
        aj = view(a, :, j)
        r[j] = dot(aj, aj)
    end
    return r
end

function _get_dists!(r::Array{Float64,2}, a::Array{Float64,2})
    m, n = size(a)
    mul!(r, a', a)
    sa2 = sumsq_percol(a)
    @inbounds for j = 1:n
        for i = 1:(j - 1)
            @inbounds r[i, j] = r[j, i]
        end
        r[j, j] = 0
        sa2j = sa2[j]
        @simd for i = (j + 1):n
            r[i, j] = sqrt(max(sa2[i] + sa2j - 2 * r[i, j], 0))
        end
    end
    r
end

function get_dists(datapoints::Array{Float64,2})
    n = size(datapoints, 2)
    r = Matrix{Float64}(undef, n, n)
    _get_dists!(r, datapoints)
    diameter = maximum(r)
    r, diameter
end

function get_neighbours(dists::Array{Float64,2}, threshold::Float64,
        num_points::Int64)
    neighbours = [BitSet() for i = 1:num_points]
    for i = 1:num_points
        @inbounds @simd for j = i+1:num_points
            @inbounds if dists[j, i] <= threshold
                push!(neighbours[i], j)
                push!(neighbours[j], i)
            end
        end
    end
    neighbours
end

################################
### CORE ALGORITHM FUNCTIONS ###
################################

function get_best_cluster(unclustered, dists, neighbours, threshold,
        candidate_clusters, candidate_diameters, calculated_dict)
    best_cluster_size = zero(Int64)
    best_diameter = zero(Float64)
    best_cluster = BitSet()
    for seed in unclustered
        if length(candidate_clusters[seed]) != 0
            candidate, candidate_diameter =
                candidate_clusters[seed], candidate_diameters[seed]
            candidate_size = length(candidate)
        else
            candidate, candidate_diameter =
                generate_candidate(seed, dists, neighbours, threshold,
                calculated_dict, candidate_clusters, candidate_diameters)
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

function generate_candidate(seed, dists, neighbours,
        threshold::Float64, calculated_dict, candidate_clusters, candidate_diameters)
    candidate_cluster = BitSet([seed])
    last_added = seed
    seed_neighbours = collect(copy(neighbours[seed]))
    diameter_cache = zeros(Float64, length(seed_neighbours))
    candidate_diameter = zero(Float64)

    while length(seed_neighbours) > 0
        if length(candidate_cluster) in SAVE_CLUSTER_VALS
            found_seed = find_cluster(candidate_cluster, calculated_dict, seed)
            if found_seed !== 0
                return candidate_clusters[found_seed], candidate_diameters[found_seed]
            end
        end

        candidate_point_index = update_diameter_cache!(diameter_cache, dists,
            seed_neighbours, last_added, threshold)
        candidate_point_value = seed_neighbours[candidate_point_index]

        if diameter_cache[candidate_point_index] < threshold
            candidate_diameter = diameter_cache[candidate_point_index]
            @inbounds push!(candidate_cluster, candidate_point_value)
            @inbounds popat!(seed_neighbours, candidate_point_index)
            @inbounds popat!(diameter_cache, candidate_point_index)
            last_added = candidate_point_value
        else
            return candidate_cluster, candidate_diameter
        end
    end

    if candidate_diameter === zero(Float64)
        return candidate_cluster, Inf
    end

    return candidate_cluster, candidate_diameter
end

function find_cluster(candidate_cluster::BitSet, calculated_dict, seed::Int64)
    calculated_index = hash(candidate_cluster)
    found_seed = get(calculated_dict, calculated_index, 0)
    if found_seed === 0
        calculated_dict[calculated_index] = seed
    end
    found_seed
end

function update_diameter_cache!(diameter_cache::Array{Float64,1},
        dists::Array{Float64,2}, seed_neighbours::Array{Int64,1},
        last_added::Int64, threshold::Float64)
    min_index = 1
    min_dist = threshold + 1
    @inbounds for (i, (cached_dist, new_dist)) in enumerate(zip(diameter_cache,
            view(dists, seed_neighbours, last_added)))
        @inbounds diameter_cache[i] = cached_dist > new_dist ? cached_dist : new_dist
        @inbounds (diameter_cache[i] < min_dist) &&
            ((min_index, min_dist) = (i, diameter_cache[i]))
    end
    min_index
end

function update_data!(unclustered::Array{Int64,1}, neighbours,candidate_clusters, best)
    setdiff!(unclustered, best)
    for unclustered_point in unclustered
        setdiff!(neighbours[unclustered_point], best)
        if !isdisjoint(candidate_clusters[unclustered_point], best)
            candidate_clusters[unclustered_point] = BitSet()
        end
    end
    nothing
end

######################
### MAIN ALGORITHM ###
######################
function QT(filename::String, _threshold::String)
    threshold, percentage, num_points = validate_input(filename, _threshold)
    names, datapoints = open(f -> read_points(f, num_points), filename)
    dists, diameter = get_dists(datapoints)

    if percentage
        threshold = threshold/100 * diameter
    elseif threshold > diameter
        @warn "Given threshold ($threshold) is greater than the diameter of the dataset ($diameter)."
    end

    neighbours = get_neighbours(dists, threshold, num_points)

    unclustered = collect(1:num_points)
    candidate_clusters = [BitSet() for i = 1:num_points]
    candidate_diameters = Vector{Float64}(undef, num_points)

    cluster_num = 1
    outfile = open("$filename.out", "w")

    while length(unclustered) > 0
        calculated_dict = Dict{UInt64,Int64}()
        best_cluster = get_best_cluster(unclustered, dists, neighbours,
            threshold, candidate_clusters, candidate_diameters, calculated_dict)

        best_cluster_array = collect(best_cluster)
        write_output(outfile, cluster_num,
            names[best_cluster_array], datapoints[:, best_cluster_array])
        update_data!(unclustered, neighbours, candidate_clusters, best_cluster)

        cluster_num += 1
    end

    close(outfile)
end


#####################
### MAIN FUNCTION ###
#####################
function main()
    arg_len = length(ARGS)
    try
        if arg_len == 2
            @time QT(ARGS[1], ARGS[2])
            exit(0)
        elseif arg_len == 1
            if ARGS[1] in ("-h", "--help")
                println(HELP_MESSAGE)
                exit(0)
            else
                println(ONE_ARG)
            end
        elseif arg_len == 0
            println(NO_ARGS)
        else
            println(TOO_MANY_ARGS)
        end
        exit(1)
    catch e
        msg = sprint(showerror, e)
        println(msg)
        exit(1)
    end
end

main()

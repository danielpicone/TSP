# This file matches all data and creates the relevant dataframe
using CSV
using DataFrames
using CPLEX
using JuMP
using LightGraphs
using SimpleWeightedGraphs
using MetaGraphs
using GraphLayout
using Combinatorics
using OffsetArrays
using GLPKMathProgInterface

include("funcs.jl")

convexity_constraint = nothing

if length(ARGS) >= 1
    println(ARGS[1])
    global n = parse(Int64,ARGS[1])
elseif !isdefined(:n)
    global n = 10
end

# solutions =  CSV.read("./data/solutions.csv")

positions = CSV.read("./data/positions_"*string(n)*"_v1.csv")
positions[:weight] = ones(n)
positions[:profit] = positions[:profit]*3

distances = get_distance(positions)

################################
# Create the RMP_τ^(I^P)
################################

# for num=10:10:40, v=1:10
#     solve_instance(num,v)
# end

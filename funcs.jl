# Hold all the function used for stabilised column generation
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

struct soln
    edge
    vertex
end

function get_distance(coord)
    n = size(coord)[1]
    D = zeros(n,n)
    for i=1:n, j=i:n
        D[i,j] = sqrt((coord[i,2]-coord[j,2])^2 + (coord[i,3] - coord[j,3])^2)
    end
    # return (D+D')./2
    return (D+D')
end

function greedy_shortest_path(D,start_node)
    current_node = start_node
    cities_visited = Set(start_node)
    path = [start_node]
    cost = 0
    for i=2:n
        new_cities = sortcols([D[current_node,:]'; collect(1:n)'])
        j = 2
        while (new_cities[2,j] in path)
            j += 1
        end
        push!(path,new_cities[2,j])
        cost += new_cities[1,j]
        current_node = Int64(new_cities[2,j])
    end
    cost += D[path[1],path[end]]

    return path, cost

end

function create_reduced_costs(rmp)
    d_hat = copy(convert(Array{Float64,2},distances))
    for i=1:n
        for j=i+1:n
            d_hat[i,j] += -new_dual[i] - new_dual[j]
            d_hat[j,i] = d_hat[i,j]
        end
    end

    return d_hat
end

function find_one_tree(d_hat, distances; num_trees = 1)
    graph = MetaGraph(n)
    for i=2:n, j=i+1:n
        add_edge!(graph, i, j)
        set_prop!(graph, i, j, :weight, d_hat[i,j])
    end
    tree = kruskal_mst(graph)
    num_v = zeros(n)
    num_v[1] = 2
    # Find number of vertices
    for edge in tree
        num_v[edge.src] += 1
        num_v[edge.dst] += 1
    end

    cost = 0; reduced_cost = 0
    for edge in tree
        if edge.src < edge.dst
            cost += distances[edge.src,edge.dst]
            reduced_cost += d_hat[edge.src,edge.dst]
        elseif edge.src > edge.dst
            cost += distances[edge.dst,edge.src]
            reduced_cost += d_hat[edge.dst,edge.src]
        else
            println("ERROR")
        end
    end
    new_edges = [d_hat[1,2:end] collect(2:n)]'
    new_edges = sortcols(new_edges)
    inds = Int64.(new_edges[2,:])
    tree_cost = []; num_v_array = []; reduced_cost = []; tree_array = []
    for i=1:num_trees
        temp_cost = cost + distances[1,inds[1]] + distances[1,inds[i+1]]
        temp_num_v = copy(num_v)
        temp_num_v[inds[1]] += 1
        temp_num_v[inds[i+1]] += 1
        temp_tree = copy(tree)
        push!(temp_tree, LightGraphs.SimpleGraphs.SimpleEdge{Int64}(1, inds[1]))
        push!(temp_tree, LightGraphs.SimpleGraphs.SimpleEdge{Int64}(1, inds[i+1]))
        push!(tree_array, temp_tree)
        temp_reduced_cost = temp_cost - sum(temp_num_v[j]*new_dual[j] for j=1:n) - γ
        if temp_reduced_cost > 0
            println("Reduced cost of a new column was positive. This was not added")
            return tree_cost, reduced_cost, tree_array, num_v_array
            break
        end
        push!(tree_cost, temp_cost)
        push!(num_v_array, temp_num_v)
        push!(reduced_cost, temp_reduced_cost)
    end

    return tree_cost, reduced_cost, tree_array, num_v_array
end

function append_new_col!(rmp_τ; objective_coefficient = nothing, constraint_coefficients = nothing)
    d_hat = create_reduced_costs(rmp_τ)
    if (objective_coefficient == nothing || constraint_coefficients == nothing)
        objective_coefficient, reduced_cost, tree_array, constraint_coefficients = find_one_tree(d_hat, distances; num_trees = 1)
    end

    if convexity_constraint == nothing
        global convexity_constraint = @constraint(rmp_τ, 0 == 1)
        push!(constraint_refs, convexity_constraint)
    end
    @variable(rmp_τ, 0 <= λ_new <= Inf, objective = objective_coefficient, inconstraints = constraint_refs, coefficients = [constraint_coefficients;1])
    setname(λ_new,string("λ[",rmp_τ.numCols - 2*n,"]"))
    append!(basic_array, 0)
    return objective_coefficient, constraint_coefficients
end

function change_objective!(rmp_τ, new_box)
    existing_columns = getobjective(rmp_τ).aff
    if length(existing_columns.vars[2*n+1:end]) >= 1
        @objective(rmp_τ, Min, sum((new_box[j]+ϵ_α)*var_dict["ub"][j] - (new_box[j]-ϵ_α)*var_dict["lb"][j] for j=1:n)
         + existing_columns.coeffs[end-num_variables:end]'*existing_columns.vars[end-num_variables:end])
    else
        @objective(rmp_τ, Min, sum((new_box[j]+ϵ_α)*var_dict["ub"][j] - (new_box[j]-ϵ_α)*var_dict["lb"][j] for j=1:n))
    end
    return new_box
end

function test_integrality(rmp)
    return mapreduce(x -> abs(getvalue(x))<0.000000001 ? true : false, &, true, Variable.(rmp, 1:(2*n)))
end
function print_solution(rmp_τ)
    # if round.(getvalue(λ),5)!=1
    #     for var in Variable.(rmp_τ, 1:rmp_τ.numCols)
    #         println(var, ": ", getvalue(var))
    #     end
    # end
    relaxed_solution_vars = Variable.(rmp_τ,1:rmp_τ.numCols)[abs.(rmp_τ.colVal) .> 1e-10]
    for sol in relaxed_solution_vars
        println(sol, " has a value of ", getvalue(sol))
    end
end
function in_box(dual, center, box)
    return dual < center-box ? center-box : (dual > center+box ? center+box : dual)
end
function get_variable_values(rmp_τ)
    return Variable.(rmp_τ,2*n+1:rmp_τ.numCols) |> getvalue
end
function check_basic!(rmp_τ, basic_array)
    for (index, solution_value) in enumerate(get_variable_values(rmp_τ))
        if round(solution_value,10)==0
            basic_array[index] += 1
        else
            basic_array[index] = 0
        end
    end
    return basic_array
end
function change_box_size(box_size, num_times_changed)
    if num_times_changed <= 1
        box_size /= 2
    elseif num_times_changed >= 8
        box_size *= 2
    end
    return box_size
end
function change_box_size!(ϵ_α, num_times_changed)
    ϵ_α = min(max(change_box_size(ϵ_α, num_times_changed), min_box_size), max_box_size)
    return ϵ_α
end
function heuristic_solution(rmp_τ)
    distance_edge_array = []; edge_array = []
    for i=1:n,j=i+1:n
        push!(distance_edge_array, (Edge(i,j),distances[i,j]))
    end
    sort!(distance_edge_array, lt = (x,y) -> x[2] > y[2], rev = true)
    relaxed_solution_vars = Variable.(rmp_τ,1:rmp_τ.numCols)[abs.(rmp_τ.colVal) .> 1e-10]
    edge_set = Set()
    for edges in relaxed_solution_vars
        [push!(edge_set,edge) for edge in columns[linearindex(edges)-2*n][1]]
    end
    for edge in edge_set
        value = 0
        for edges in relaxed_solution_vars
            if edge in columns[linearindex(edges)-2*n][1]
                value += edges |> getvalue
            end
        end
        push!(edge_array, (edge,value))
    end
    sort!(edge_array, lt = (x,y) -> x[2] > y[2])
    graph = MetaGraph(n)
    index = 1; solution_cost = 0
    while (!(all(length.(inneighbors.(graph,collect(1:n))).==2)) && index <= length(edge_array))
        if (length(inneighbors(graph,edge_array[index][1].src))<2 && length(inneighbors(graph,edge_array[index][1].dst))<2)
            add_edge!(graph, edge_array[index][1])
            solution_cost += distances[edge_array[index][1].src, edge_array[index][1].dst]
        end
        index += 1
    end
    index = 1
    if ne(graph) < n
        while !(all(length.(inneighbors.(graph,collect(1:n))).==2))
            if index > length(distance_edge_array)
                println("Heuristic failed to find a solution")
                return false
            end
            if (length(inneighbors(graph,distance_edge_array[index][1].src))<2 && length(inneighbors(graph,distance_edge_array[index][1].dst))<2)
                add_edge!(graph, distance_edge_array[index][1])
                solution_cost += distance_edge_array[index][2]
            end
            index += 1
        end
    end
    return edge_set, edge_array, graph, solution_cost
end


function solve_instance(num, version; ϵ = 0.0001, exclude_columns = false,
    increase_exclude_bound = false,
    change_box_size_option = false,
    use_heuristic = false,
    include_initial_columns = 0,
    num_since_basic = 20)
    global n = num
    positions = CSV.read("./data/positions_"*string(n)*"_v"*string(version)*".csv")
    global distances = get_distance(positions)
    global τ = 1
    println("Create the RMP")
    α_tilde = zeros(n)
    global ϵ_α = 10
    min_box_size = 10; max_box_size = 100000

    # Get cost of inital tour
    initial_cost = sum(distances[i,i+1] for i=1:n-1) + distances[1,n]
    global columns = Dict()

    global rmp_τ = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 0, CPXPARAM_Preprocessing_Dual=-1, CPXPARAM_MIP_Display = 2))

    global var_dict = Dict()
    var_dict["ub"] = @variable(rmp_τ, η_ub[1:n]>=0)
    var_dict["lb"] = @variable(rmp_τ, η_lb[1:n]>=0)

    # Add the objective function
    @objective(rmp_τ, Min, sum((α_tilde[i] + ϵ_α)*η_ub[i] - (α_tilde[i] - ϵ_α)*η_lb[i] for i=1:n))

    @constraint(rmp_τ, vertex_constraints[i=1:n], η_ub[i] - η_lb[i] == 2)
    global constraint_refs = Vector{ConstraintRef}(0)
    for i=1:n
        push!(constraint_refs, vertex_constraints[i])
    end

    lower_bound = []; upper_bound = []; objective_array = []
    duals = []; global basic_array = []; global constraint_array = []; global non_basic_array = []
    println("Solving the CGP now")

    # Variables
    num_times_v_lb_changed = 0
    v_ub = Inf; v_lb = -Inf; gap = Inf
    solved = false
    global convexity_constraint = nothing
    temp = 0
    push!(duals, α_tilde)


    global new_dual = α_tilde

    global num_variables = 0
    start_time = time()
    for i=1:include_initial_columns
        new_values = greedy_shortest_path(distances,i)
        append_new_col!(rmp_τ, objective_coefficient = new_values[2], constraint_coefficients = 2*ones(n))
        num_variables += 1
    end
    while ((gap > ϵ) & !solved)
        solve(rmp_τ)
        objective_value = getobjectivevalue(rmp_τ)
        if test_integrality(rmp_τ)
            v_ub = objective_value
        end
        append!(upper_bound, v_ub)
        new_dual = getdual(vertex_constraints)
        if convexity_constraint != nothing
            global γ = getdual(convexity_constraint)
        else
            global γ = 0
        end
        if (exclude_columns & (τ % 100 == 0))
            check_basic!(rmp_τ, basic_array)
            for var in Variable.(rmp_τ, 2*n+1:rmp_τ.numCols)[basic_array .> num_since_basic]
                if increase_exclude_bound
                    if MathProgBase.getconstrmatrix(rmp_τ |> internalmodel)[:, var.col] in constraint_array
                        num_since_basic += 5
                    end
                end
                println("THIS RUNS")
                # setlowerbound(var, 0.0)
                setupperbound(var, 0.0)
            end
        end
        if τ == 100
            for var in Variable.(rmp_τ, 2*n+1:2*n+1+include_initial_columns)
                setupperbound(var, 0.0)
            end
        end
        # TODO: fix up change_objective!, it is allowing duals to move too much
        new_box = in_box.(new_dual, duals[τ], ϵ_α)
        new_dual = change_objective!(rmp_τ, new_box)
        push!(duals, new_dual)

        d_hat = create_reduced_costs(rmp_τ)
        new_objective_coefficients, reduced_cost, tree_array, new_constraint_coefficients = find_one_tree(d_hat, distances; num_trees = 1)
        if (use_heuristic && (test_integrality(rmp_τ) && (τ % 100 == 0)))
            heur_solution = heuristic_solution(rmp_τ)
            if heur_solution != false
                append_new_col!(rmp_τ, objective_coefficient = heur_solution[4], constraint_coefficients = 2*ones(n))
            end
            push!(non_basic_array, rmp_τ.numCols)
        elseif !isempty(new_objective_coefficients)
            for (index,new_obs) in enumerate(new_objective_coefficients)
                objective_coefficient, constraint_coefficients  = append_new_col!(rmp_τ; objective_coefficient = new_obs, constraint_coefficients = new_constraint_coefficients[index])
                if all(constraint_coefficients.==2)
                    println("An integral solution was found with cost: ", objective_coefficient)
                    solved = true
                    break
                end
            end
        end
        # Remove heuristic solutions from basic solutions
        for variable in non_basic_array
            setupperbound(Variable.(rmp_τ, variable),0.0)
        end
        temp = copy(reduced_cost)
        if length(reduced_cost) == 0
            break
        end
        append!(lower_bound, v_lb)
        push!(constraint_array, MathProgBase.getconstrmatrix(rmp_τ |> internalmodel)[:, rmp_τ.numCols])
        if τ>1
            v_lb = max(v_lb, objective_value + reduced_cost[1])
            gap = (v_ub-v_lb)/abs(v_lb)
        else
            v_lb = objective_value + reduced_cost[1]
        end
        if ((τ > 100) & change_box_size_option)
            num_times_v_lb_changed = lower_bound[end-9:end] |> unique |> length
            ϵ_α = change_box_size!(ϵ_α, num_times_v_lb_changed)
        end
        println("At iteration ", τ, " gap is: ",gap)
        for (index,tree) in enumerate(tree_array)
            columns[rmp_τ.numCols - 2*n] = (tree, new_objective_coefficients[index])
        end
        τ+=1; num_variables += 1
        # if τ > 12
        #     break
        # end
    end

    solve(rmp_τ)
    new_dual = getdual(vertex_constraints)
    push!(duals, new_dual)
    end_time = time()
    println("Completed ",τ," iterations in: ", end_time - start_time, " seconds")
    # relaxed_solution_vars = Variable.(rmp_τ,1:rmp_τ.numCols)[abs.(rmp_τ.colVal) .> 1e-10]
    # for sol in relaxed_solution_vars
    #     println(sol, " has a value of ", getvalue(sol))
    # end

    return duals, lower_bound, upper_bound

    integer_solution = heuristic_solution(rmp_τ)
end

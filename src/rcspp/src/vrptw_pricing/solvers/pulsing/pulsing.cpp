#include "rcspp/vrptw_pricing/solver/pulsing/pulsing.h"
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <tuple>

namespace rcspp::vrptw_pricing::pulsing{

using rcspp::graph::adjacency_matrix::VertexId;;

namespace {

struct ResourceState{
    Cost cost;
    Time time;
    Demand demand;
};

struct PulsingState{

    struct ResultComparator{
        bool operator()(const Result& a, const Result& b){
            return a.cost < b.cost; 
        }
    };

    std::vector<VertexId> path;
    std::vector<bool> is_visited;

    std::vector<Time> depart_times;
    std::vector<Cost> depart_costs;
  
    Cost min_primal_cost = INF_COST;
    std::vector<Result> results;
    ResultComparator result_comp;

    bool bounding_mode = false;
    Cost min_cost_bound = INF_COST;

    PulsingState(VertexId vertices_num, VertexId source, bool bounding_mode = false)
    : path{source}, is_visited(vertices_num), depart_times(vertices_num), depart_costs(vertices_num), bounding_mode(bounding_mode)
    {
        is_visited[source] = true;
    }
};

struct Workspace{

    struct CostBound {
        Time time;
        Cost cost_bound;
    };

    using Neighbors = std::vector<VertexId>;
    using CostBounds = std::vector<CostBound>;
    
    std::vector<Neighbors> sorted_successors_per_vertex;
    std::vector<CostBounds> sorted_bounds_per_vertex;

    SolveParams solve_params{};

    Workspace(VertexId vertex_count, const SolveParams& solve_params = {})
    : sorted_successors_per_vertex(vertex_count), sorted_bounds_per_vertex(vertex_count), solve_params(solve_params)
    {}
};

/**
 * @brief Implements the core recursive backtracking search of the pulsing algorithm.
 * This function explores one path at a time in a depth-first manner. Upon arriving at a
 * new vertex, it performs a series of pruning checks:
 * 
 * 1.  Feasibility Pruning: Checks for violations of time windows and vehicle capacity.
 * 
 * 2.  Target Check: If the target is reached, the path is a candidate solution, and the
 *     primal bound (best solution cost) is updated.
 * 
 * 3.  Rollback Pruning: Checks if a "shortcut" path from the grand-predecessor is better.
 * 
 * 4.  Cost-Bound Pruning: Uses pre-computed bounds to prune paths that cannot possibly improve
 *     upon the current best solution.
 *
 * If the path is not pruned, the function recursively calls itself for each valid successor.
 * After the recursive calls return, it backtracks by removing the current vertex from the path.
 *
 * @param context The problem context.
 * @param workspace The workspace containing pre-computed data like sorted successors.
 * @param pulsing_state The current state of the search (e.g., current path, visited nodes).
 * @param arrival_state The resource state (cost, time, demand) upon arrival at the current vertex.
 */
void pulsing(const Context& context, const Workspace& workspace, PulsingState& pulsing_state, const ResourceState& arrival_state){

    const Graph& graph = context.graph;
    
    if(pulsing_state.path.empty()){
        throw std::logic_error("Path shoud not be empty.");
    }

    VertexId current_vertex_id = pulsing_state.path.back();
    const Vertex& cur_vertex = graph.get_vertex(current_vertex_id);

    // --- Feasibity pruning ---

    // Whether violate time window
    if(arrival_state.time > cur_vertex.due_time){
        return;
    }

    Time effective_arrival_time = std::max(arrival_state.time, cur_vertex.ready_time);

    // Whether overload
    if(arrival_state.demand + cur_vertex.demand > context.capacity){
        return;
    }

    // --- Update primal and return if reach target ---
    if(current_vertex_id == context.target){
        Cost candidate_primal_cost = arrival_state.cost + cur_vertex.cost;
        
        if(!pulsing_state.bounding_mode && candidate_primal_cost < context.cost_threshold){
            std::vector<Result>& results = pulsing_state.results;
            results.push_back({candidate_primal_cost, pulsing_state.path});
            std::push_heap(results.begin(), results.end(), pulsing_state.result_comp);
            while ((results.size() > context.max_paths))
            {
                std::pop_heap(results.begin(), results.end(), pulsing_state.result_comp);
                results.pop_back();
            }
        }
        
        pulsing_state.min_primal_cost = std::min(pulsing_state.min_primal_cost, candidate_primal_cost);
        
        if(pulsing_state.bounding_mode){
            pulsing_state.min_cost_bound = std::min(pulsing_state.min_cost_bound, candidate_primal_cost);
        }

        // The reason why check "size >= 2" is that
        // the source and target can be the same, and we shouldn't stop pulsing
        // if the target is the fist vertex.
        if(pulsing_state.path.size() >= 2){
            return;
        }
    }

    // Rollback pruning
    // Current path is start -> ... -> u -> v -> w (current vertex), 
    // is path start -> ... -> u -> w not worse than it?
    int64_t path_size = pulsing_state.path.size();
    if(path_size >= 3){
        VertexId u = pulsing_state.path[path_size - 3];
        const Edge* edge_ptr = graph.try_get_edge(u, current_vertex_id);
        if(edge_ptr != nullptr){
            Cost candidate_arrival_cost = pulsing_state.depart_costs[u] + edge_ptr->cost;
            Time candidate_effective_arrival_time = std::max(cur_vertex.ready_time, pulsing_state.depart_times[u] + edge_ptr->time);
            if(candidate_arrival_cost <= arrival_state.cost
                && candidate_effective_arrival_time <= effective_arrival_time)
            {
                return;
            }
        }
    }

    // --- Cost bound pruning ---

    // Retrieve cost bound given cumuative time
    const Workspace::CostBounds& sorted_bounds = workspace.sorted_bounds_per_vertex[current_vertex_id];
    auto it = std::lower_bound(sorted_bounds.begin(), sorted_bounds.end(), effective_arrival_time,
                    [](const Workspace::CostBound& bound, Time effective_arrival_time){
                        return bound.time > effective_arrival_time;
                    }
                ); // Find first bound whose time <= arrival time
    Cost cost_bound = -INF_COST;
    if(it != sorted_bounds.end()){
        cost_bound = arrival_state.cost + it->cost_bound;

        if(pulsing_state.bounding_mode && !workspace.solve_params.exact_bounding){
            pulsing_state.min_cost_bound = std::min(pulsing_state.min_cost_bound, cost_bound);
            return;
        }
    }

    // Primal cost pruning
    if(cost_bound >= pulsing_state.min_primal_cost){
        return;
    }

    // Threshold pruning (e.g. reduced cost constant item)
    if(cost_bound >= context.cost_threshold && !pulsing_state.bounding_mode){
        return;
    }

    // --- Update pulse state ---
    Cost depart_cost = arrival_state.cost + cur_vertex.cost;
    Time depart_time = effective_arrival_time + cur_vertex.service_time;
    Demand depart_demand = arrival_state.demand + cur_vertex.demand;

    pulsing_state.depart_costs[current_vertex_id] = depart_cost;
    pulsing_state.depart_times[current_vertex_id] = depart_time;

    // --- Expand to successors ---
    for(VertexId u: workspace.sorted_successors_per_vertex[current_vertex_id]){
        if(pulsing_state.is_visited[u] && u != context.target){
            // The reason why check "next vertex != target" is that
            // the source and target can be the same.
            // So we should expand to target even if itself is visited before.
            continue;
        }

        const Edge& edge_props = graph.get_edge(current_vertex_id, u);

        pulsing_state.path.push_back(u);
        pulsing_state.is_visited[u] = true;

        ResourceState next_arrival_state{depart_cost + edge_props.cost, depart_time + edge_props.time, depart_demand};
        pulsing(context, workspace, pulsing_state, next_arrival_state);

        pulsing_state.path.pop_back();
        pulsing_state.is_visited[u] = false;
    }
}

/**
 * @brief Pre-computes heuristic lower bounds on path costs for pruning.
 * This function improves the effectiveness of the pulsing algorithm's pruning. It discretizes
 * the time resource into segments. For each time segment and for each vertex 'u', it runs the
 * `pulsing` algorithm in a special "bounding mode" to find the minimum possible cost to travel from 'u'
 * to the target. These pre-computed costs are stored in the workspace and are later used by the
 * main `pulsing` function for more aggressive cost-bound pruning.
 * @param t_upper The upper bound of the time horizon for bounding.
 * @param t_lower The lower bound of the time horizon.
 * @param t_delta The time discretization step.
 * @param context The problem context.
 * @param workspace The workspace where the computed bounds will be stored.
 */
void bounding(Time t_upper, Time t_lower, Time t_delta, const Context& context, Workspace& workspace){
    const Graph& graph = context.graph;

    Time t = t_upper - t_delta;
    while (t > t_lower + 1){
        for(const Vertex& u: graph.vertices()){
            if(u.id == context.source || u.id == context.target){
                continue;
            }
            // Time window check
            if(t > u.due_time){
                workspace.sorted_bounds_per_vertex[u.id].push_back({t, INF_COST});
                continue;
            }
            if(t + t_delta <= u.ready_time){
                continue;
            }

            bool bounding_mode = true;
            PulsingState pulsing_state(graph.vertices_num(), u.id, bounding_mode);
            ResourceState arrival_state{0.0, t, 0.0};
            pulsing(context, workspace, pulsing_state, arrival_state);

            workspace.sorted_bounds_per_vertex[u.id].push_back({t, pulsing_state.min_cost_bound});
        }
        t -= t_delta;
    }
    
}

/**
 * @brief Pre-sorts the successors for each vertex based on a heuristic traversal cost.
 * For each vertex 'u', this function iterates through its successors 'v' and calculates a
 * simple heuristic cost (edge cost + successor vertex cost). It then sorts the successors in
 * ascending order of this cost. This heuristic guides the depth-first search of the pulsing
 * algorithm to explore more promising paths first, which can lead to the earlier discovery of
 * good solutions and thus more effective cost-bound pruning.
 * @param context The problem context, containing the graph.
 * @param workspace The workspace object where the sorted successors will be stored.
 */
void sort_successors(const Context& context, Workspace& workspace){

    struct SortableVertex{
        // bool is_target;
        Cost traversal_cost;
        VertexId vertex_id;
    };
    
    const Graph& graph = context.graph;

    for(const Vertex& u: graph.vertices()){

        std::vector<SortableVertex> successors;
        successors.reserve(graph.out_degree(u));

        for(const Vertex& v: graph.successors(u)){
            Cost traversal_cost = graph.get_edge(u, v).cost + v.cost;
            successors.push_back({
                traversal_cost,
                v.id
            });
        }

        std::sort(successors.begin(), successors.end(), 
            [](const SortableVertex& a, const SortableVertex& b){
                return a.traversal_cost < b.traversal_cost;
        });


        Workspace::Neighbors& sorted_successors = workspace.sorted_successors_per_vertex[u.id];
        sorted_successors.reserve(successors.size());
        for(const SortableVertex& succ: successors){
            sorted_successors.push_back(succ.vertex_id);
        }
    }
}

}

/**
 * @brief The main public entry point to run the VRPTW pricing pulsing algorithm.
 * This function orchestrates the entire solution process. It first sets up a workspace,
 * then performs pre-computation steps by calling `sort_successors` and, optionally, the `bounding`
 * procedure. It then initializes the main search state and kicks off the recursive `pulsing`
 * process from the source vertex. Finally, it collects, sorts, and returns all feasible paths
 * found that are better than the specified cost threshold.
 * @param context The problem context, defining the graph, source, target, and constraints.
 * @param solve_params Parameters to control the algorithm's behavior, such as bounding settings.
 * @return A `std::vector<Result>` containing the best paths found by the algorithm, sorted by cost.
 */
std::vector<Result> solve(const Context& context, const SolveParams& solve_params){
    const VertexId source = context.source;
    const VertexId target = context.target;
    const Graph& graph = context.graph;

    if(!graph.is_vertex_valid(source) || !graph.is_vertex_valid(target)){
        return {};
    }

    Workspace workspace(graph.vertices_num(), solve_params);
    sort_successors(context, workspace);

    // --- Bounding ---
    Time t_upper = graph.get_vertex(target).due_time;
    Time t_lower = graph.get_vertex(source).ready_time;
    Time t_delta = solve_params.bounding_t_delta;
    if(t_upper > t_lower && t_delta > 0){
        double t_seg_num = std::ceil((t_upper - t_lower) / t_delta);
        if(t_seg_num >= 2){
            t_delta = std::ceil((t_upper - t_lower) / t_seg_num); // Adjust delta to make length of segments even
            bounding(t_upper, t_lower, t_delta, context, workspace);
        }
    }

    // --- Solve initial state ---

    PulsingState pulsing_state(graph.vertices_num(), source);
    ResourceState arrival_state{0.0, 0.0, 0.0};
    pulsing(context, workspace, pulsing_state, arrival_state);
    std::sort_heap(pulsing_state.results.begin(), pulsing_state.results.end(), pulsing_state.result_comp);
    return pulsing_state.results;
}

}
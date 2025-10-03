#include "rcspp/vrptw_pricing/solver/labeling/labeling.h"
#include <unordered_set>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rcspp::vrptw_pricing::labeling{

using rcspp::graph::adjacency_matrix::VertexId;;

namespace{

using LabelId = int64_t;
inline const LabelId INVALID_LABEL_ID = -1;
inline const VertexId INVALID_VERTEX_ID = -1;

struct ResourceState{
    Cost cost;
    Time time;
    Demand demand;
};

struct Label{
public:
    LabelId id = INVALID_LABEL_ID;
    LabelId predecessor_id = INVALID_LABEL_ID;
    VertexId vertex_id = INVALID_VERTEX_ID;
    int64_t depth = 0;

    Cost arrival_cost = 0;
    Time arrival_time = 0;
    Demand arrival_demand = 0;

    Cost depart_cost = 0;
    Time depart_time = 0;
    Demand depart_demand = 0;

    std::vector<VertexId> ng_memorized_vertices;
    bool is_dominated = false;
    
private:
    Label() = default;
    friend class LabelBuilder;
};

class LabelBuilder {
public:
    LabelBuilder(LabelId predecessor_id, VertexId vertex_id, int64_t depth) {
        instance.predecessor_id = predecessor_id;
        instance.vertex_id = vertex_id;
        instance.depth = depth;
    }

    LabelBuilder& with_arrival(Cost c, Time t, Demand d) {
        instance.arrival_cost = c;
        instance.arrival_time = t;
        instance.arrival_demand = d;
        return *this;
    }

    LabelBuilder& with_depart(Cost c, Time t, Demand d) {
        instance.depart_cost = c;
        instance.depart_time = t;
        instance.depart_demand = d;
        return *this;
    }

    LabelBuilder& with_ng_memory(const std::vector<VertexId>& ng_memorized_vertices){
        instance.ng_memorized_vertices = ng_memorized_vertices;
        return *this;
    }
    
    Label build() {
        return instance;
    }

private:
    Label instance;
};

class LabelContainer{
public:
    const Label& get_label(LabelId label_id) const{
        if(!is_label_valid(label_id)){
            throw std::logic_error("Label does not exist.");
        }
        return labels_[label_id];
    }

    Label& get_label(LabelId label_id){
        return const_cast<Label&>( 
            static_cast<const LabelContainer*>(this)->get_label(label_id)
        );
    }

    void add_label(Label& label){
        label.id = labels_.size();
        labels_.push_back(label);
    }


private:
    std::vector<Label> labels_;

    bool is_label_valid(LabelId label_id) const{
        return label_id >= 0 && label_id < labels_.size();
    }
};

struct Workspace{

    struct CostBound {
        Time time;
        Cost cost_bound;
    };

    using Neighbors = std::vector<VertexId>;
    using CostBounds = std::vector<CostBound>;
    
    std::vector<Neighbors> ng_neighbors_per_vertex;
    std::vector<CostBounds> sorted_bounds_per_vertex;

    SolveParams solve_params{};

    Workspace(VertexId vertices_num, const SolveParams& solve_params = {})
    : ng_neighbors_per_vertex(vertices_num), sorted_bounds_per_vertex(vertices_num), solve_params(solve_params)
    {}
};

enum struct Dominance{
    Dominating,
    Dominated,
    NonDominance,
    Equal
};

struct LabelingState{

    using LabelIdPool = std::unordered_set<LabelId>;

    struct SortableLabel{
        Cost f_cost;
        int64_t depth;
        LabelId label_id;

        bool operator>(const SortableLabel& other) const{
            if(f_cost != other.f_cost){
                return f_cost > other.f_cost;
            }
            if(depth != other.depth){
                return depth < other.depth;
            }
            return label_id > other.label_id;
        }
    };

    struct ResultComparator {
        bool operator()(const Result& a, const Result& b) const {
            return a.cost < b.cost;
        }
    };

    LabelContainer labels;
    std::vector<LabelIdPool> pareto_labels_per_vertex;
    std::priority_queue<SortableLabel, std::vector<SortableLabel>, std::greater<SortableLabel>> label_queue;

    Cost min_primal_cost = INF_COST;
    std::vector<Result> results;
    ResultComparator result_comp;

    bool bounding_mode = false;
    Cost min_cost_bound = INF_COST;

    LabelingState(VertexId vertices_num, bool bounding_mode = false)
        : pareto_labels_per_vertex(vertices_num), bounding_mode(bounding_mode) {}
};

/**
 * @brief Pre-computes the NG-Relaxation neighbors for each vertex in the graph.
 * For each vertex 'u', this function identifies the 'k' successors that are cheapest to travel to,
 * based on the sum of edge cost and successor vertex cost. These pre-computed neighbors are
 * stored in the Workspace and used later for the partial elementarity check (NG-Relaxation pruning),
 * which helps to reduce cycles in paths. The number of neighbors 'k' is determined by
 * `solve_params.ng_neighbor_size`.
 * @param context The problem context, containing the graph.
 * @param workspace The workspace object where the computed NG neighbors will be stored.
 */
void init_ng_neighbors(const Context& context, Workspace& workspace){
    if(workspace.solve_params.ng_neighbor_size <= 0){
        return;
    }

    struct SortableVertex{
        Cost traversal_cost;
        VertexId vertex_id;
    };

    const Graph& graph = context.graph;

    for(const Vertex& u: graph.vertices()){
        std::vector<SortableVertex> neighbors;
        neighbors.reserve(graph.out_degree(u));

        for(const Vertex& v: graph.successors(u)){
            Cost traversal_cost = graph.get_edge(u, v).cost + v.cost;
            neighbors.push_back({traversal_cost, v.id});
        }

        std::sort(neighbors.begin(), neighbors.end(), 
            [](const SortableVertex& a, const SortableVertex& b){
                return a.traversal_cost < b.traversal_cost;
            }
        );

        Workspace::Neighbors& ng_neighbors = workspace.ng_neighbors_per_vertex[u.id];
        int64_t ng_size = std::max<int64_t>(0, std::min<int64_t>(workspace.solve_params.ng_neighbor_size, graph.out_degree(u)));
        ng_neighbors.reserve(ng_size);
        for(int64_t i = 0; i < neighbors.size() && i < ng_size; ++i){
            ng_neighbors.push_back(neighbors[i].vertex_id);
        }
        std::sort(ng_neighbors.begin(), ng_neighbors.end());
    }
}

/**
 * @brief Reconstructs the sequence of visited vertices from a given final label.
 * This function backtracks from the provided `label` through its chain of predecessors
 * using the `predecessor_id`. It collects the vertex ID at each step, rebuilding the
 * path that led to this label.
 * @param labeling_state The current state of the algorithm, used to access the global label container.
 * @param label The final label from which to reconstruct the path.
 * @return A `std::vector<VertexId>` representing the path from the source vertex to the target vertex.
 */
std::vector<VertexId> construct_path(const LabelingState& labeling_state, const Label& label){
    std::vector<VertexId> path;
    LabelId cur_label_id = label.id;
    while(true){
        const Label& cur_label = cur_label_id == label.id ? label: labeling_state.labels.get_label(cur_label_id);
        path.push_back(cur_label.vertex_id);
        if(cur_label.depth <= 0){
            break;
        }
        cur_label_id = cur_label.predecessor_id;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

/**
 * @brief Compares two labels to determine their dominance relationship.
 * Label 'a' dominates label 'b' if it
 * arrives at the same vertex with resources that are at least as good in every dimension (cost,
 * effective arrival time, demand, and NG-memorized set) and strictly better in at least one
 * dimension. The NG-set comparison is based on set inclusion.
 * @param a The first label for comparison.
 * @param b The second label for comparison.
 * @param context The problem context, needed for vertex-specific data like ready times.
 * @return A `Dominance` enum indicating whether 'a' dominates 'b', is dominated by 'b',
 * is non-dominant, or is equal.
 */
Dominance check_dominance(const Label& a, const Label& b, const Context& context){
    if(a.vertex_id != b.vertex_id){
        return Dominance::NonDominance;
    }

    bool has_better_dim = false;
    bool has_worse_dim = false;

    // Compare cost
    if(a.arrival_cost < b.arrival_cost){
        has_better_dim = true;
    }else if(a.arrival_cost > b.arrival_cost){
        has_worse_dim = true;
    }

    // Compare effective arrival time
    Time ready_time = context.graph.get_vertex(a.vertex_id).ready_time;
    Time a_effective_arrival_time = std::max(a.arrival_time, ready_time);
    Time b_effective_arrival_time = std::max(b.arrival_time, ready_time);

    if(a_effective_arrival_time < b_effective_arrival_time){
        has_better_dim = true;
    }else if(a_effective_arrival_time > b_effective_arrival_time){
        has_worse_dim = true;
    }

    // Compare demand
    if(a.arrival_demand < b.arrival_demand){
        has_better_dim = true;
    }else if(a.arrival_demand > b.arrival_demand){
        has_worse_dim = true;
    }

    // Compare ng-memorized vertices
    if(std::includes(b.ng_memorized_vertices.begin(), b.ng_memorized_vertices.end(),
                    a.ng_memorized_vertices.begin(), a.ng_memorized_vertices.end()))
    {
        if(a.ng_memorized_vertices.size() < b.ng_memorized_vertices.size()){
            has_better_dim = true;
        }
    }else if(std::includes(a.ng_memorized_vertices.begin(), a.ng_memorized_vertices.end(),
                            b.ng_memorized_vertices.begin(), b.ng_memorized_vertices.end()))
    {
        if(a.ng_memorized_vertices.size() > b.ng_memorized_vertices.size()){
            has_worse_dim = true;
        }
    }else{
        has_better_dim = true;
        has_worse_dim = true;
    }

    if(has_better_dim && has_worse_dim){
        return Dominance::NonDominance;
    }else if(has_better_dim && !has_worse_dim){
        return Dominance::Dominating;
    }else if(!has_better_dim && has_worse_dim){
        return Dominance::Dominated;
    }else{
        return Dominance::Equal;
    }
}

/**
 * @brief Expands a path from a predecessor label to a new successor vertex.
 * This is the main workhorse of the labeling algorithm. It performs the following steps:
 * 
 * 1.  Calculates the new resource state (cost, time, demand) upon arriving at successor `vertex`.
 * 
 * 2.  Performs a series of pruning checks:
 *     
 *     - Feasibility: Violates time windows or capacity.
 *     
 *     - NG-Relaxation: Prunes if `vertex` is in the predecessor's NG-path memory.
 *     
 *     - Rollback Pruning: Checks if a direct path from the grand-predecessor is better.
 * 
 * 3.  If a path to the target is found, updates the primal bound (best solution cost).
 * 
 * 4.  Calculates a heuristic f-cost (A* search) and prunes if it exceeds the current best solution.
 * 
 * 5.  Constructs a candidate label and checks it for dominance against existing Pareto-optimal labels.
 * 
 * 6.  If the candidate label is not dominated, it is added to the priority queue and the set of
 *     Pareto labels for the given vertex.
 * @param predecessor The label representing the path so far.
 * @param vertex The successor vertex to expand to.
 * @param context The problem context.
 * @param workspace The workspace containing pre-computed data.
 * @param labeling_state The current algorithm state, which is modified by adding new labels.
 */
void expand_label(const Label& predecessor, const Vertex& vertex, const Context& context, const Workspace& workspace, LabelingState& labeling_state){

    const Graph& graph = context.graph;
    LabelContainer& labels = labeling_state.labels;

    // Partial elementarity pruning using ng-relaxation
    bool is_ng_memorized = std::binary_search(predecessor.ng_memorized_vertices.begin(), predecessor.ng_memorized_vertices.end(), vertex.id);
    if(is_ng_memorized && vertex.id != context.target){
        // The reason why check "vertex != target" is that
        // the source and target can be the same.
        // So we should expand to target even if itself is visited before (and thus may already be memorized).
        return;
    }

    Cost arrival_cost = predecessor.depart_cost + (predecessor.vertex_id == vertex.id ? 0: graph.get_edge(predecessor.vertex_id, vertex.id).cost);
    Time arrival_time = predecessor.depart_time + (predecessor.vertex_id == vertex.id ? 0: graph.get_edge(predecessor.vertex_id, vertex.id).time);
    Time effective_arrival_time = std::max(arrival_time, vertex.ready_time);
    Demand arrival_demand = predecessor.depart_demand;

    Cost depart_cost = arrival_cost + vertex.cost;
    Time depart_time = effective_arrival_time + vertex.service_time;
    Demand depart_demand = arrival_demand + vertex.demand;

    // --- Feasibility pruning ---
    if(arrival_time > vertex.due_time){
        return;
    }
    if(arrival_demand + vertex.demand > context.capacity){
        return;
    }

    // --- Rollback pruning ---
    // If path is start -> ... -> u -> v (predecessor) -> w (vertex), 
    // is path start -> ... -> u -> w (vertex) not worse than it?
    if(predecessor.depth >= 1){
        const Label& label_u = labels.get_label(predecessor.predecessor_id);
        const Edge* edge_ptr = graph.try_get_edge(label_u.vertex_id, vertex.id);
        if(edge_ptr != nullptr){
            Cost candidate_arrival_cost = label_u.depart_cost + edge_ptr->cost;
            Time candidate_effective_arrival_time = std::max(label_u.depart_time + edge_ptr->time, vertex.ready_time);
            if(candidate_arrival_cost <= arrival_cost
                && candidate_effective_arrival_time <= effective_arrival_time)
            {
                return;
            }
        }
    }

    // --- Construct candidate label ---
    Label candidate_label = LabelBuilder(predecessor.id, vertex.id, predecessor.depth + 1)
                                .with_arrival(arrival_cost, arrival_time, arrival_demand)
                                .with_depart(depart_cost, depart_time, depart_demand)
                                .build();

    // --- Update primal and return if reach target ---
    // We should NOT update primal cost only when the top label of queue reach target.
    // The case indicates procedure found optimal and will terminate, so primal cost is not used at all. 
    if(vertex.id == context.target){
        Cost candidate_primal_cost = depart_cost;

        if(!labeling_state.bounding_mode && candidate_primal_cost < context.cost_threshold){
            // Collect path
            std::vector<Result>& results = labeling_state.results;
            results.push_back({
                candidate_primal_cost,
                construct_path(labeling_state, candidate_label) 
            });
            std::push_heap(results.begin(), results.end(), labeling_state.result_comp);
            while(results.size() > context.max_paths){
                std::pop_heap(results.begin(), results.end(), labeling_state.result_comp);
                results.pop_back();
            }
        }

        labeling_state.min_primal_cost = std::min(labeling_state.min_primal_cost, candidate_primal_cost);

        if(labeling_state.bounding_mode){
            labeling_state.min_cost_bound = std::min(labeling_state.min_cost_bound, candidate_primal_cost);
        }

        // The reason why check "depth >= 1" is that
        // the source and target can be the same, and we shouldn't stop pulsing
        // if the target is the fist vertex.
        if(candidate_label.depth >= 1){
            return;
        }
    }

    // --- F-cost pruning ---

    // Calculate f-cost
    Cost f_cost = -INF_COST;
    if(vertex.id == context.target && candidate_label.depth >= 1){
        // The reason why check "depth >= 1" is that
        // the source and target can be the same, and we shouldn't just use depart cost
        // as f-cost if the target is the fist vertex (so depth of is 0).
        f_cost = depart_cost;
    }else{
        const Workspace::CostBounds& sorted_bounds = workspace.sorted_bounds_per_vertex[vertex.id];
        auto it = std::lower_bound(sorted_bounds.begin(), sorted_bounds.end(), effective_arrival_time,
                        [](const Workspace::CostBound& bound, Time effective_arrival_time){
                            return bound.time > effective_arrival_time;
                        }
                    ); // Find first bound whose time <= arrival time
        if(it != sorted_bounds.end()){
            f_cost = arrival_cost + it->cost_bound;

            if(labeling_state.bounding_mode && !workspace.solve_params.exact_bounding){
                labeling_state.min_cost_bound = std::min(labeling_state.min_cost_bound, f_cost);
                return;
            }
        }
    }

    // f-cost pruning
    if(f_cost >= labeling_state.min_primal_cost || (f_cost >= context.cost_threshold && !labeling_state.bounding_mode)){
        return;
    }


    // --- Construct ng-memorized vertices of candidate label ---

    // Intersection with predecessor
    const Workspace::Neighbors& ng_neighbors = workspace.ng_neighbors_per_vertex[vertex.id];
    std::vector<VertexId>& ng_memorized_vertices = candidate_label.ng_memorized_vertices;
    ng_memorized_vertices.reserve(ng_neighbors.size());
    std::set_intersection(
        predecessor.ng_memorized_vertices.begin(), predecessor.ng_memorized_vertices.end(),
        ng_neighbors.begin(), ng_neighbors.end(),
        std::back_inserter(ng_memorized_vertices)
    );

    // Union with current vertex
    auto it = std::lower_bound(ng_memorized_vertices.begin(), ng_memorized_vertices.end(), vertex.id);
    if(it == ng_memorized_vertices.end() || *it != vertex.id){
        ng_memorized_vertices.insert(it, vertex.id);
    }

    // --- Dominance pruning ---
    bool dominated = false;
    bool equal = false;
    auto& pareto_label_ids = labeling_state.pareto_labels_per_vertex[vertex.id];
    for(auto it=pareto_label_ids.begin(); it != pareto_label_ids.end();){
        Label& active_label = labels.get_label(*it);
        Dominance dominance = check_dominance(candidate_label, active_label, context);
        switch (dominance)
        {
            case Dominance::Dominating:
                it = pareto_label_ids.erase(it);
                active_label.is_dominated = true;
                break;
            case Dominance::NonDominance:
                ++it;
                break;
            case Dominance::Dominated:
                dominated = true;
                break;
            case Dominance::Equal:
                equal = true;
                break;
            default:
                break;
        }
        if(dominated || equal){
            break;
        }
    }

    if(dominated || equal){
        return;
    }

    // --- Add candidate label to container, pareto label pool of associated vertex and label queue ---
    
    labeling_state.labels.add_label(candidate_label);
    pareto_label_ids.insert(candidate_label.id);
    labeling_state.label_queue.push({f_cost, candidate_label.depth, candidate_label.id});
}

/**
 * @brief Creates and prepares the initial label at the source vertex to start the algorithm.
 * This function bootstraps the labeling process. It avoids duplicating logic by
 * constructing an artificial predecessor label at the source vertex and then calling `expand_label`
 * on the source vertex itself. This allows the initialization to reuse all the feasibility checking,
 * f-cost calculation, and queuing logic contained within `expand_label`.
 * @param source_vertex The source vertex of the graph.
 * @param context The problem context.
 * @param workspace The workspace.
 * @param labeling_state The algorithm state to be initialized.
 * @param arrival_state The initial resource state at the source (typically all zeros).
 */
void prepare_init_label(const Vertex& source_vertex, const Context& context, const Workspace& workspace, LabelingState& labeling_state, const ResourceState& arrival_state){
    int64_t depth = -1;
    Label predecessor = LabelBuilder(INVALID_LABEL_ID, source_vertex.id, depth)
                        .with_arrival(arrival_state.cost, arrival_state.time, arrival_state.demand)
                        .with_depart(arrival_state.cost, arrival_state.time, arrival_state.demand)
                        .build();
    expand_label(predecessor, source_vertex, context, workspace, labeling_state);
}

/**
 * @brief Executes the main loop of the A*-based labeling algorithm.
 * This function continuously extracts the most promising label (one with the lowest f-cost)
 * from a priority queue. For each extracted label, it performs validity checks (e.g., has it been
 * dominated since it was queued?) and termination checks. If the label is valid and promising,
 * it iterates through all successor vertices and calls `expand_label` to extend the path. The loop
 * continues until the queue is empty or the termination condition (f-cost > primal bound) is met.
 * @param context The problem context.
 * @param workspace The workspace.
 * @param labeling_state The current algorithm state.
 */
void labeling(const Context& context, Workspace& workspace, LabelingState& labeling_state){

    const Graph& graph = context.graph;

    auto& labels = labeling_state.labels;
    auto& label_queue = labeling_state.label_queue;

    while (!label_queue.empty()){
        const auto top_label = label_queue.top();
        label_queue.pop();
        
        LabelId label_id = top_label.label_id;
        Cost f_cost = top_label.f_cost;
        const Label label = labels.get_label(label_id);
        const Vertex& cur_vertex = graph.get_vertex(label.vertex_id);

        // --- Dominance skip check ---
        if(label.is_dominated){
            continue;
        }
        
        // --- F-cost termination check ---
        if(f_cost >= labeling_state.min_primal_cost
            || (f_cost >= context.cost_threshold && !labeling_state.bounding_mode)){
            return;
        }

        // --- Expand to successors ---
        for(const Vertex& u: graph.successors(cur_vertex)){
            expand_label(label, u, context, workspace, labeling_state);
        }
    }
}

/**
 * @brief Pre-computes heuristic lower bounds on path costs for the f-cost calculation.
 * This function improves the A* heuristic by calculating tighter lower bounds. It discretizes
 * the time resource into segments. For each time segment and for each vertex 'u', it runs a separate
 * labeling procedure from 'u' to the target to find the minimum possible cost. These
 * pre-computed costs are stored in the workspace and later used by `expand_label` to get a more
 * accurate f-cost, leading to more aggressive pruning.
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
            LabelingState labeling_state(graph.vertices_num(), bounding_mode);
            const ResourceState arrival_state{0, t, 0};
            prepare_init_label(u, context, workspace, labeling_state, arrival_state);
            labeling(context, workspace, labeling_state);

            workspace.sorted_bounds_per_vertex[u.id].push_back({t, labeling_state.min_cost_bound});
        }
        t -= t_delta;
    }
}

}

/**
 * @brief The main public entry point to run the VRPTW pricing labeling algorithm.
 * This function orchestrates the entire solution process. It sets up a workspace,
 * initializes NG-neighbors, and optionally runs the `bounding` procedure to generate
 * A* heuristic cost bounds. It then initializes the main labeling state, prepares the initial label
 * at the source vertex, and executes the `labeling` algorithm. Finally, it collects and returns
 * all feasible paths found that are better than the specified cost threshold.
 * @param context The problem context, defining the graph, source, target, and constraints.
 * @param solve_params Parameters to control the algorithm's behavior, such as NG-neighbor size.
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
    init_ng_neighbors(context, workspace);

    // --- Bounding ---
    Time t_upper = graph.get_vertex(source).due_time;
    Time t_lower = graph.get_vertex(target).ready_time;
    Time t_delta = solve_params.bounding_t_delta;
    if(t_upper > t_lower && t_delta > 0){
        double t_seg_num = std::ceil((t_upper - t_lower) / t_delta);
        if(t_seg_num >= 2){
            t_delta = std::ceil((t_upper - t_lower) / t_seg_num); // Adjust delta to make length of segments even
            bounding(t_upper, t_lower, t_delta, context, workspace);
        }
    }

    // --- Solve initial state ---
    
    LabelingState labeling_state(graph.vertices_num());
    // Prepare initial label
    const ResourceState arrival_state{0, 0, 0};
    prepare_init_label(context.graph.get_vertex(context.source), context, workspace, labeling_state, arrival_state);
    labeling(context, workspace, labeling_state);
    std::sort_heap(labeling_state.results.begin(), labeling_state.results.end(), labeling_state.result_comp);
    return labeling_state.results;
}


}
#pragma once

#include "rcspp/graph/types.h"
#include "rcspp/graph/adjacency_matrix.h"
#include <limits>
#include <cstdint>
#include <vector>

namespace rcspp::vrptw_pricing
{

using rcspp::graph::adjacency_matrix::VertexId;

using Cost = double;
using Time = double;
using Demand = double;

inline const Cost INF_COST = std::numeric_limits<Cost>::infinity();

struct Vertex{
    VertexId id;
    Time ready_time;
    Time due_time;
    Time service_time;
    Demand demand;
    Cost cost;
};

struct Edge{
    Cost cost;
    Time time;
};

using Graph = rcspp::graph::adjacency_matrix::Graph<Vertex, Edge>;
    
struct Context{
    Graph graph;

    VertexId source;
    VertexId target;

    Demand capacity;
    Cost cost_threshold = INF_COST;
    int64_t max_paths = std::numeric_limits<int64_t>::max();

    Context(VertexId vertices_num)
        : graph(Graph::create(vertices_num))
    {}

    Context(VertexId vertices_num, VertexId source, VertexId target, Demand capacity)
        : graph(Graph::create(vertices_num)), source(source), target(target), capacity(capacity)
    {}
};

struct Result{
    Cost cost;
    std::vector<VertexId> path;
};


}

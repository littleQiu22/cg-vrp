#include "rcspp/graph/types.h"
#include "rcspp/common/iterator_range.h"
#include <iterator>
#include <optional>
#include <vector>
#include <string>

namespace rcspp::graph::adjacency_matrix{

using rcspp::graph::Degree;
using rcspp::graph::VertexId;

// Forward declaration
template <typename Vertex, typename Edge>
class Graph;

namespace detail{

template<typename Vertex, typename Edge>
class SuccessorIterator{
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Vertex;
    using difference_type = std::ptrdiff_t;
    using pointer = const Vertex*;
    using reference = const Vertex&;

    SuccessorIterator(const Graph<Vertex, Edge>* graph, VertexId u, VertexId v)
        : graph_(graph), source_(u), current_target_(v){}

    reference operator*() const{
        return graph_->get_vertex(current_target_);
    }

    SuccessorIterator& operator++(){
        for(VertexId u = current_target_ + 1; u < graph_->vertices_num(); ++u){
            if(graph_->try_get_edge(source_, u)){
                current_target_ = u;
                return *this;
            }
        }
        current_target_ = graph_->vertices_num();
        return *this;
    }

    SuccessorIterator operator++(int) {
        SuccessorIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    friend bool operator==(const SuccessorIterator& a, const SuccessorIterator& b){
        return a.source_ == b.source_ && a.current_target_ == b.current_target_;
    }

    friend bool operator!=(const SuccessorIterator& a, const SuccessorIterator& b){
        return !(a == b);
    }

private:
    const Graph<Vertex, Edge>* graph_;
    VertexId source_;
    VertexId current_target_;
};

template<typename Vertex, typename Edge>
class PredecessorIterator{
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Vertex;
    using difference_type = std::ptrdiff_t;
    using pointer = const Vertex*;
    using reference = const Vertex&;

    PredecessorIterator(const Graph<Vertex, Edge>* graph, VertexId u, VertexId v)
        : graph_(graph), current_source_(u), target_(v){}

    reference operator*() const{
        return graph_->get_vertex(current_source_);
    }

    PredecessorIterator& operator++(){
        for(VertexId u = current_source_ + 1; u < graph_->vertices_num(); ++u){
            if(graph_->try_get_edge(u, target_)){
                current_source_ = u;
                return *this;
            }
        }
        current_source_ = graph_->vertices_num();
        return *this;
    }

    PredecessorIterator operator++(int) {
        PredecessorIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    friend bool operator==(const PredecessorIterator& a, const PredecessorIterator& b){
        return a.current_source_ == b.current_source_ && a.target_ == b.target_;
    }

    friend bool operator!=(const PredecessorIterator& a, const PredecessorIterator& b){
        return !(a == b);
    }

private:
    const Graph<Vertex, Edge>* graph_;
    VertexId current_source_;
    VertexId target_;
};

template<typename Vertex, typename Edge>
class VertexIterator {
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type        = Vertex;
    using difference_type   = std::ptrdiff_t;
    using pointer           = const Vertex*;
    using reference         = const Vertex&;

    explicit VertexIterator(const Graph<Vertex, Edge>* graph, VertexId u = 0) 
        : graph_(graph), current_vertex_id_(u) {}

    reference operator*() const {
        return graph_->get_vertex(current_vertex_id_);
    }

    VertexIterator& operator++() {
        ++current_vertex_id_;
        return *this;
    }
    
    VertexIterator operator++(int) {
        VertexIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    friend bool operator==(const VertexIterator& a, const VertexIterator& b) {
        return a.current_vertex_id_ == b.current_vertex_id_;
    }

    friend bool operator!=(const VertexIterator& a, const VertexIterator& b) {
        return !(a == b);
    }

private:
    const Graph<Vertex, Edge>* graph_;
    VertexId current_vertex_id_;
};

}

/**
 * @brief A graph representation using adjacency matrix as data container.
 * @tparam Vertex The type of data associated with each vertex.
 * @tparam Edge The type of data associated with each edge.
 *
 */
template<typename Vertex, typename Edge>
class Graph{
public:
    using SuccessorIterator = detail::SuccessorIterator<Vertex, Edge>;
    using SuccessorRange = IteratorRange<SuccessorIterator>;
    using PredecessorIterator = detail::PredecessorIterator<Vertex, Edge>;
    using PredecessorRange = IteratorRange<PredecessorIterator>;
    using VertexIterator = detail::VertexIterator<Vertex, Edge>;
    using VertexRange = IteratorRange<VertexIterator>;

    static Graph create(VertexId vertices_num) {
        if (vertices_num < 0) {
            throw std::invalid_argument("Vertex count cannot be negative.");
        }
        return Graph(vertices_num);
    }

    VertexId vertices_num() const{
        return vertices_num_;
    }

    // --- Edge operations ---

    bool set_edge(VertexId source, VertexId target, const Edge& edge)
    {
        if(!is_vertex_valid(source) || !is_vertex_valid(target) || source == target){
            return false;
        }

        matrix_[source][target] = edge;
        ++out_degrees_[source];
        ++in_degrees_[target];
        return true;
    }

    bool remove_edge(VertexId source, VertexId target){
        if(!is_vertex_valid(source) || !is_vertex_valid(target)){
            return false;
        }

        matrix_[source][target] = std::nullopt;
        --out_degrees_[source];
        --in_degrees_[target];
        return true;
    }

    bool remove_edge(const Vertex& source, const Vertex& target){
        return remove_edge(source.id, target.id);
    }

    const Edge* try_get_edge(VertexId source, VertexId target) const{
        if(!is_vertex_valid(source) || !is_vertex_valid(target)){
            return nullptr;
        }
        const auto& edge_opt = matrix_[source][target];
        if (edge_opt) {
            return &(*edge_opt);
        }
        return nullptr;
    }

    const Edge* try_get_edge(const Vertex& source, const Vertex& target) const{
        return try_get_edge(source.id, target.id);
    }

    const Edge& get_edge(VertexId source, VertexId target) const{
        const Edge* edge_ptr = try_get_edge(source, target);
        if(edge_ptr == nullptr){
            throw std::logic_error("Edge (" + std::to_string(source) + ", " + std::to_string(target) + ") does not exists.");
        }
        return *edge_ptr;
    }

    const Edge& get_edge(const Vertex& source, const Vertex& target) const{
        return get_edge(source.id, target.id);
    }

    // --- Vertex operations ---

    bool set_vertex(const Vertex& vertex){
        VertexId vertex_id = vertex.id;
        if(!is_vertex_valid(vertex.id)){
            return false;
        }
        vertices_[vertex_id] = vertex;
        return true;
    }

    const Vertex* try_get_vertex(VertexId vertex_id) const{
        if(!is_vertex_valid(vertex_id)){
            return nullptr;
        }
        return &vertices_[vertex_id];
    }

    Vertex* try_get_vertex(VertexId vertex_id){
        return const_cast<Vertex*>(static_cast<const Graph<Vertex, Edge>*>(this)->try_get_vertex(vertex_id));
    }

    const Vertex& get_vertex(VertexId vertex_id) const{
        const Vertex* vertex_ptr = try_get_vertex(vertex_id);
        if(vertex_ptr == nullptr){
            throw std::logic_error("Vertex " + std::to_string(vertex_id) + " does not exist.");
        }
        return *vertex_ptr;
    }

    Vertex& get_vertex(VertexId vertex_id){
        return const_cast<Vertex&>(static_cast<const Graph<Vertex, Edge>*>(this)->get_vertex(vertex_id));
    }

    SuccessorRange successors(VertexId vertex_id) const{
        VertexId first_target = vertices_num();
        if(is_vertex_valid(vertex_id)){
            for(VertexId u = 0; u < vertices_num(); ++u){
                if(matrix_[vertex_id][u]){
                    first_target = u;
                    break;
                }
            }
        }

        return SuccessorRange(
            SuccessorIterator(this, vertex_id, first_target),
            SuccessorIterator(this, vertex_id, vertices_num())
        );
    }

    SuccessorRange successors(const Vertex& vertex) const{
        return successors(vertex.id);
    }

    PredecessorRange predecessors(VertexId vertex_id) const{
        VertexId first_source = vertices_num();
        if(is_vertex_valid(vertex_id)){
            for(VertexId u = 0; u < vertices_num(); ++u){
                if(matrix_[u][vertex_id]){
                    first_source = u;
                    break;
                }
            }
        }

        return PredecessorRange(
            PredecessorIterator(this, first_source, vertex_id),
            PredecessorIterator(this, vertices_num(), vertex_id)
        );
    }

    PredecessorRange predecessors(const Vertex& vertex) const{
        return predecessors(vertex.id);
    }

    VertexRange vertices() const {
        return VertexRange(
            VertexIterator(this, 0),
            VertexIterator(this, vertices_num())
        );
    }

    bool is_vertex_valid(VertexId vertex_id) const {
        return vertex_id >= 0 && vertex_id < vertices_num_;
    }

    bool is_vertex_valid(const Vertex& vertex) const {
        return is_vertex_valid(vertex.id);
    }

    Degree out_degree(VertexId vertex_id) const{
        if(!is_vertex_valid(vertex_id)){
            return 0;
        }
        return out_degrees_[vertex_id];
    }

    Degree in_degree(VertexId vertex_id) const{
        if(!is_vertex_valid(vertex_id)){
            return 0;
        }
        return in_degrees_[vertex_id];
    }

    Degree out_degree(const Vertex& vertex) const{
        return out_degrees_[vertex.id];
    }

    Degree in_degree(const Vertex& vertex) const{
        return in_degrees_[vertex.id];
    }

private:

    explicit Graph(VertexId vertices_num)
        :vertices_num_(vertices_num),
        matrix_(vertices_num, std::vector<std::optional<Edge>>(vertices_num)),
        vertices_(vertices_num), in_degrees_(vertices_num), out_degrees_(vertices_num)
    {}

    VertexId vertices_num_;
    std::vector<std::vector<std::optional<Edge>>> matrix_;
    std::vector<Vertex> vertices_;
    std::vector<Degree> out_degrees_;
    std::vector<Degree> in_degrees_;

};

}
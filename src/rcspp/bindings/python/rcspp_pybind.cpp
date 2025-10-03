#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/optional.h>
#include "rcspp/vrptw_pricing/problem.h"
#include "rcspp/vrptw_pricing/solver/labeling/labeling.h"
#include "rcspp/vrptw_pricing/solver/pulsing/pulsing.h"

namespace nb = nanobind;

void bind_vrptw_pricing(nb::module_& m) {
    
    using namespace rcspp::vrptw_pricing;

    nb::class_<Vertex>(m, "Vertex")
        .def(nb::init<VertexId, Time, Time, Time, Demand, Cost>(), nb::arg("id"), nb::arg("ready_time"), nb::arg("due_time"), nb::arg("service_time"), nb::arg("demand"), nb::arg("cost"));

    nb::class_<Edge>(m, "Edge")
        .def(nb::init<Cost, Time>(), nb::arg("cost"), nb::arg("time"));

    nb::class_<Graph>(m, "Graph")
        .def("set_edge", &Graph::set_edge, nb::arg("source"), nb::arg("target"), nb::arg("edge"))
        .def("set_vertex", &Graph::set_vertex, nb::arg("vertex"))
        .def("get_edge", nb::overload_cast<VertexId, VertexId>(&Graph::get_edge, nb::const_), nb::arg("source"), nb::arg("target"), nb::rv_policy::copy)
        .def("get_vertex", nb::overload_cast<VertexId>(&Graph::get_vertex, nb::const_), nb::arg("vertex_id"), nb::rv_policy::copy);


    nb::class_<Context>(m, "Context")
        .def(nb::init<VertexId>(), nb::arg("vertices_num"))
        .def(nb::init<VertexId, VertexId, VertexId, Demand>(), nb::arg("vertices_num"), nb::arg("source"), nb::arg("target"), nb::arg("capacity"))
        .def_ro("graph", &Context::graph)
        .def_rw("capacity", &Context::capacity)
        .def_rw("cost_threshold", &Context::cost_threshold)
        .def_rw("max_paths", &Context::max_paths);

    nb::class_<Result>(m, "Result")
        .def_ro("cost", &Result::cost)
        .def_ro("path", &Result::path);

    nb::class_<labeling::SolveParams>(m, "LabelingParams")
        .def(nb::init<int64_t, bool, Time>(), nb::arg("ng_neighbor_size"), nb::arg("exact_bounding"), nb::arg("bounding_t_delta"));

    nb::class_<pulsing::SolveParams>(m, "PulsingParams")
        .def(nb::init<bool, double>(), nb::arg("exact_bounding"), nb::arg("bounding_t_delta"));

    m.def("labeling", &labeling::solve, nb::arg("context"), nb::arg("params"), nb::rv_policy::move);
    m.def("pulsing", &pulsing::solve, nb::arg("context"), nb::arg("params"), nb::rv_policy::move);

}

NB_MODULE(_rcspp, m){
    m.doc() = "A C++ library for Resource-Constrained-Shortest-Path-Problem.";
    
    nb::module_ vrptw_module = m.def_submodule("vrptw_pricing", "Solvers for Vehicle Routing with Time Windows pricing.");
    bind_vrptw_pricing(vrptw_module);
}
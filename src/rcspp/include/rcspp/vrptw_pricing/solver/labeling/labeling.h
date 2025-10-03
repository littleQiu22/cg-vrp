#pragma once

#include "rcspp/vrptw_pricing/problem.h"
#include <cstdint>

namespace rcspp::vrptw_pricing::labeling{

struct SolveParams{
    int64_t ng_neighbor_size = 8;
    bool exact_bounding = true;
    Time bounding_t_delta = 10;
};

std::vector<Result> solve(const Context& context, const SolveParams& solve_params = SolveParams{});

}
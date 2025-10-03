#pragma once

#include "rcspp/graph/types.h"
#include "rcspp/vrptw_pricing/problem.h"

namespace rcspp::vrptw_pricing::pulsing{

struct SolveParams{
    bool exact_bounding = true;
    double bounding_t_delta = 10.0;
};

std::vector<Result> solve(const Context& context, const SolveParams& solve_params = SolveParams{});

}
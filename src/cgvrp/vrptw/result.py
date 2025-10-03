from cgvrp.vrptw.route import Route

from dataclasses import dataclass
from typing import Sequence, Tuple

lp_value = float

@dataclass
class Result:
    """
    Stores the final result of the Vehicle Routing Problem with Time Window solution.

    This is a data class that encapsulates the final solution, including the cost,
    the dual bound, the optimality gap, and the routes that make up the solution.

    Attributes:
        cost (float): The final objective value of the solution (primal value).
        bound (float): The final dual bound calculated during the column generation process.
        gap (float): The final optimality gap, calculated as |cost - bound| / |cost|.
        routes (Sequence[Tuple[Route, lp_value]]): A sequence of tuples, where each tuple
            contains a Route object and its corresponding value (lp_value) in the
            final linear programming solution.
    """

    cost: float
    bound: float
    gap: float

    routes: Sequence[Tuple[Route, lp_value]]
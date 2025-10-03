from cgvrp.vrptw.problem import Problem, vertex_id
from cgvrp.vrptw.result import Result
from cgvrp.vrptw.route import Route, RouteSource
from cgvrp.vrptw.model import Model
from cgvrp import vrptw_pricing

from typing import Set, List, Tuple, Optional, Dict, Literal
from dataclasses import dataclass
import time

@dataclass
class PricedRoute:
    """
    A data class to store a route found by the pricing subproblem.

    Contains both the route itself and its calculated reduced cost.

    Attributes:
        reduced_cost (float): The reduced cost of the route.
        route (Route): The route object.
    """

    reduced_cost: float
    route: Route

class ColumnGeneration:
    """
    Main class that uses the column generation algorithm for the VRPTW.

    This class manages the main loop of the algorithm, alternating between solving
    the Restricted Master Problem (RMP) and the pricing subproblem to generate
    new columns (routes) until a convergence criterion is met.
    """

    def __init__(self, problem: Problem) -> None:
        self._validate_problem(problem)
        self.problem = problem

        self.model_: Optional[Model] = None
        self.pricing_ctx: Optional[vrptw_pricing.Context] = None
        self.vertex_to_ctx_idx: Optional[Dict[vertex_id, int]] = None
        self.ctx_idx_to_vertex: Optional[Dict[int, vertex_id]] = None
        
        self.iter: Optional[int] = None
        self.primal: Optional[float] = None
        self.bound: Optional[float] = None
        self.gap: Optional[float] = None
        self.lp_time = 0.0
        self.pricing_time = 0.0

    def _validate_problem(self, problem: Problem):
        """
        Performs basic validations on the provided problem instance.

        Raises:
            ValueError: If a depot has not been set in the problem.
        """

        if problem.depot is None:
            raise ValueError("Depot is not set.")
        
    def _print_log_header(self):
        header = f"{'Iter':<6} {'Status':<20} {'Primal':<12} {'Bound':<12} {'Gap(%)':<12} {'LP_Time(s)':<12} {'Pricing_Time(s)':<14} {'#Columns':<12}"
        print(header)
        print("-" * len(header))

    def _print_iter_log(self, status: Literal["Solving LP", "Pricing", "Iteration Complete"], end: str = "\n"):
        primal_str = f"{self.primal:.2f}" if status != "Solving LP" else ""
        bound_str = f"{self.bound:.2f}" if status == "Iteration Complete" else ""
        gap_str = f"{self.gap*100:.2f}" if status == "Iteration Complete" else ""
        lp_time_str = f"{self.lp_time:.2f}" if status != "Solving LP" else ""
        pricing_time_str = f"{self.pricing_time:.2f}" if status == "Iteration Complete" else ""
        columns_str = f"{len(self.model.route_to_var)}" if status == "Iteration Complete" else ""
        
        log_line = f"{self.iter:<6} {status:<20} {primal_str:<12} {bound_str:<12} {gap_str:<12} {lp_time_str:<12} {pricing_time_str:<14} {columns_str:<12}"
        print(log_line, end=end, flush=True)

    def solve(self, gap_limit: float=1e-2, pricing_method: Literal["labeling", "pulsing"]="pulsing")->Result:
        """
        Executes the main column generation loop.

        The method iteratively solves the master problem's LP relaxation, uses the dual prices
        to solve the pricing subproblem, adds new routes (columns) with negative reduced
        cost to the master problem, and recalculates the lower bound. The loop terminates
        when the gap between the primal and dual bounds is within the specified limit.

        Args:
            gap_limit (float, optional): The target optimality gap to terminate. Defaults to 1e-2.
            pricing_method (Literal["labeling", "pulsing"], optional): The algorithm to use for the
                pricing subproblem. Defaults to "pulsing".

        Returns:
            Result: A Result object containing the final solution.
        """

        self.model = self._create_model()
        self.pricing_ctx, self.vertex_to_ctx_idx, self.ctx_idx_to_vertex = self._create_pricing_context()

        self.iter = 0
        self.primal: Optional[float] = None
        self.bound = 0.0
        self.gap: Optional[float] = None
        self.result: Optional[Result] = None
        self.lp_time = 0.0
        self.pricing_time = 0.0

        # Calculate time delta when pricing algorithm runs bounding procedure
        total_edge_time = 0.0
        total_service_time = 0.0
        edge_num = self.problem.graph.number_of_edges()
        vertex_num = self.problem.graph.number_of_nodes()
        bounding_t_delta = 0.0
        for _, _, data in self.problem.graph.edges(data=True):
            total_edge_time += data["time"]
        for u, data in self.problem.graph.nodes(data=True):
            total_service_time += data["service_time"]
        if edge_num > 0:
            bounding_t_delta += total_edge_time / edge_num
        if vertex_num > 0:
            bounding_t_delta += total_service_time / vertex_num

        self._print_log_header()
        # Main loop
        while True:
            self.iter += 1

            # Solve Restricted-Master-Problem
            _st = time.time()
            self._print_iter_log("Solving LP", end="\r")
            self.model.solve()
            self.primal = self.model.obj
            self.lp_time += time.time() - _st

            # Pricing
            _st = time.time()
            self._print_iter_log("Pricing", end="\r")
            priced_routes = self.pricing(pricing_method, bounding_t_delta)
            self.pricing_time += time.time() - _st

            # Estimate lower bound and check termination condition
            min_reduced_cost = 0.0
            for priced_route in priced_routes:
                min_reduced_cost = min(min_reduced_cost, priced_route.reduced_cost)
                self.model.add_route(priced_route.route)
            self.bound = max(self.bound, self.primal + min_reduced_cost * self.problem.vehicle_num)
            self.gap = abs(self.primal - self.bound) / max(1e-10, abs(self.primal))
            self._print_iter_log("Iteration Complete")

            if self.gap <= gap_limit:
                routes = []
                for route, lp_value in self.model.route_to_value.items():
                    if lp_value > 0:
                        routes.append((route, lp_value))
                self.result = Result(
                    cost=self.primal,
                    bound=self.bound,
                    gap=self.gap,
                    routes=routes
                )
                print("\n=== Reach gap limit ===")
                return self.result

    def pricing(self, pricing_method: Literal["labeling", "pulsing"]="pulsing", bounding_t_delta: float=0.0)->List[PricedRoute]:
        """
        Solves the pricing subproblem to find routes with negative reduced costs.

        This method updates the costs of the pricing graph using the dual variables from the
        master problem. It then calls the specified C++ pricing algorithm (labeling or pulsing)
        to find the shortest path with resource constraints.
        Pulsing algorithm is guaranteed to return elementary path, while labeling algorithm may return non-elementary path.

        Args:
            pricing_method (Literal["labeling", "pulsing"], optional): The algorithm to use. Defaults to "pulsing".
            bounding_t_delta (float, optional): A time delta used for bounding in the pricing algorithm. Defaults to 0.0 (no bounding).

        Returns:
            List[PricedRoute]: A list of new routes with negative reduced costs.
        """

        # Update reduced cost constant and vertex cost according to dual
        for u, constr in self.model.vertex_cover_constrs.items():
            ctx_idx = self.vertex_to_ctx_idx[u]
            vertex_data = self.problem.graph.nodes[u]
            self.pricing_ctx.graph.set_vertex(vrptw_pricing.Vertex(
                id=ctx_idx,
                ready_time=vertex_data["ready_time"],
                due_time=vertex_data["due_time"],
                service_time=vertex_data["service_time"],
                demand=vertex_data["demand"],
                cost=-constr.pi
            ))

        reduced_cost_constant = -self.model.vehicle_limit_constr.pi
        self.pricing_ctx.cost_threshold = -reduced_cost_constant

        # Run pricing algorithm
        if pricing_method == "pulsing":
            routes = vrptw_pricing.pulsing(context=self.pricing_ctx, params=vrptw_pricing.PulsingParams(
                exact_bounding=True, bounding_t_delta=bounding_t_delta
            ))
        else:
            routes = vrptw_pricing.labeling(context=self.pricing_ctx, params=vrptw_pricing.LabelingParams(
                exact_bounding=True, bounding_t_delta=bounding_t_delta, ng_neighbor_size=1
            ))

        pricing_routes: List[PricedRoute] = []
        for route in routes:
            reduced_cost = route.cost + reduced_cost_constant
            if reduced_cost >= 0:
                continue
            actual_cost = self._cal_cost(route.path)
            pricing_route = PricedRoute(
                reduced_cost=reduced_cost,
                route=Route(RouteSource.PRICER, cost=actual_cost, path=route.path)
            )
            pricing_routes.append(pricing_route)
        return pricing_routes
        

    def _create_model(self):
        """
        Creates and initializes the Restricted Master Problem (RMP).

        Initializes the set covering and vehicle limit constraints and populates the model
        with an initial set of routes to ensure feasibility.

        Returns:
            Model: The initialized RMP model.
        """

        model = Model()
        model.init_vehicle_limit_constr(self.problem.vehicle_num)
        for u in self.problem.graph.nodes:
            if u == self.problem.depot:
                continue
            model.init_cover_constr(u)

        routes = self._gen_init_routes()
        for route in routes:
            model.add_route(route)
        return model

    def _create_pricing_context(self):
        """
        Creates the context required by the C++ pricing subproblem solver.

        This involves mapping the problem's graph vertices and edges to the data structures
        expected by the C++ `vrptw_pricing` extension, including setting up capacities,
        time windows, and demands.

        Returns:
            Tuple[vrptw_pricing.Context, Dict[vertex_id, int], Dict[int, vertex_id]]: A tuple containing
                the pricing context, a mapping from problem vertex IDs to context indices, and the
                reverse mapping.
        """

        vertex_to_ctx_idx: Dict[vertex_id, int] = {}
        ctx_idx_to_vertex: Dict[int, vertex_id] = {}
        for u in self.problem.graph.nodes:
            ctx_idx = len(vertex_to_ctx_idx)
            vertex_to_ctx_idx[u] = ctx_idx
            ctx_idx_to_vertex[ctx_idx] = u

        context = vrptw_pricing.Context(
            vertices_num=self.problem.graph.number_of_nodes(),
            source=vertex_to_ctx_idx[self.problem.depot],
            target=vertex_to_ctx_idx[self.problem.depot],
            capacity=self.problem.vehicle_capacity
        )
        context.max_paths = 50

        for (u, v, data) in self.problem.graph.edges(data=True):
            u_ctx_idx = vertex_to_ctx_idx[u]
            v_ctx_idx = vertex_to_ctx_idx[v]
            context.graph.set_edge(u_ctx_idx, v_ctx_idx, vrptw_pricing.Edge(
                cost=data["cost"],
                time=data["time"]
            ))
        for u, vertex_data in self.problem.graph.nodes(data=True):
            ctx_idx = vertex_to_ctx_idx[u]
            context.graph.set_vertex(vrptw_pricing.Vertex(
                id=ctx_idx,
                ready_time=vertex_data["ready_time"],
                due_time=vertex_data["due_time"],
                service_time=vertex_data["service_time"],
                demand=vertex_data["demand"],
                cost=0.0
            ))

        return (context, vertex_to_ctx_idx, ctx_idx_to_vertex)

    def _gen_init_routes(self):
        """
        Generate initial routes for Restricted-Master-Problem (RMP). 
        
        Purpose of initial routes:
        1. [Necessary] help column generation (CG) succeed:
            - If Mater Problem (MP) is feasible, initial routes should make RMP always feasible, otherwise CG will fail
            - Initial routes shouldn't change the fact that OPTIMAL_VALUE(MP) = OPTIMAL_VALUE(RMP)
        2. [Optional] make column generation procedure more efficient:
            - try find routes that appears in the optimal solution of VRP  

        For purpose #1, we can construct artificial routes that:
        - have high objective cost (e.g. big-M)
        - doesn't increment the amount of used vehicles

        For purpose #2, we can use lightweight constructive heuritiscs:
        - Clarke and Wright Savings Method
        """

        routes: List[Route] = []
        depot = self.problem.depot

        # --- Artificial routes ---

        # Calculate big-M
        bigM = 0
        for (_, _, data) in self.problem.graph.edges(data=True):
            bigM += data["cost"]
        
        # Route with big-M cost
        for u in self.problem.graph.nodes:
            if u == depot:
                continue
            routes.append(Route(
                source=RouteSource.ARTIFICIAL,
                path=[depot, u, depot],
                cost=bigM
            ))
        
        # --- Routes by Clarke and Wright savings heuristic ---
        for route in self._cw_savings():
            routes.append(route)
        
        return routes
        
    def _cw_savings(self):
        """
        Implements the Clarke and Wright savings heuristic.

        This method generates a set of initial feasible routes by iteratively merging routes
        that result in the greatest cost savings.

        Returns:
            List[Route]: A list of routes generated by the heuristic.
        """

        graph = self.problem.graph
        depot = self.problem.depot

        routes: Set[Route] = set()
        vertex_to_route: Dict[vertex_id, Route] = {}

        # Initial routes
        for u in graph.nodes:
            if u == depot:
                continue
            path = [depot, u, depot]
            if not self._is_path_valid(path):
                continue
            route = Route(
                RouteSource.CW_SAVINGS,
                0.0, # Re-calculated it after all routes are determined 
                path
            )
            routes.add(route)
            vertex_to_route[u] = route

        # Calculate savings of merging routes
        savings: list[Tuple[vertex_id, vertex_id, float]] = [] # element: (vertex 1, vertex 2, saving if from vertex 1 to vertex 2)
        for u, v, edge_data in graph.edges(data=True):
            if (u, depot) not in graph.edges or (depot, v) not in graph.edges:
                continue
            saving = graph.edges[u, depot]["cost"] + graph.edges[depot, v]["cost"] - edge_data["cost"]
            if saving <= 0:
                continue
            savings.append((u, v, saving))
        
        # Merging routes
        savings.sort(key=lambda x: x[2], reverse=True)
        for u, v, saving in savings:
            if u not in vertex_to_route or v not in vertex_to_route:
                continue

            route_u = vertex_to_route[u]
            route_v = vertex_to_route[v]

            if route_u == route_v:
                continue

            if route_u.path[-2] != u or route_v.path[1] != v:
                continue

            merged_path = route_u.path[:-1] + route_v.path[1:]
            if not self._is_path_valid(merged_path):
                continue

            merged_route = Route(
                source=RouteSource.CW_SAVINGS,
                cost=0.0,
                path=merged_path
            )

            routes.remove(route_u)
            routes.remove(route_v)
            routes.add(merged_route)

            for k in merged_path:
                if k == 0:
                    continue
                vertex_to_route[k] = merged_route
        
        for route in routes:
            route.cost = self._cal_cost(route.path)
        return list(routes)

    def _is_path_valid(self, path: List[vertex_id])->bool:
        """
        Checks if a given path is feasible.

        A path is feasible if it respects vehicle capacity and all customer time windows.

        Args:
            path (List[vertex_id]): A list of vertex IDs representing the path.

        Returns:
            bool: True if the path is feasible, False otherwise.
        """

        graph = self.problem.graph

        arrival_time = 0
        cumu_demand = 0
        prev: Optional[vertex_id] = None
        for cur in path:
            if prev is not None and (prev, cur) not in graph.edges:
                return False
            
            vertex_data = graph.nodes[cur]
            if prev is not None:
                edge_data = graph.edges[prev, cur]
                arrival_time += edge_data["time"]
            if arrival_time > vertex_data["due_time"]:
                return False
            
            cumu_demand += vertex_data["demand"]
            if cumu_demand > self.problem.vehicle_capacity:
                return False
            
            arrival_time = max(vertex_data["ready_time"], arrival_time) + vertex_data["service_time"]
            prev = cur
        return True

    def _cal_cost(self, seq: List[vertex_id])->float:
        """
        Calculates the total travel cost of a sequence of vertices.

        Args:
            seq (List[vertex_id]): The sequence of vertices in the path.

        Returns:
            float: The total cost of the path.
        """

        cost = 0
        prev: Optional[vertex_id] = None
        for u in seq:
            if prev is not None:
                cost += self.problem.graph.edges[prev, u]["cost"]
            prev = u
        return cost

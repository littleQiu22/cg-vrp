from cgvrp.vrptw.problem import vertex_id
from cgvrp.vrptw.route import Route, RouteSource

import pulp

from typing import Dict, Optional

class Model:
    """
    Manages the Restricted Master Problem (RMP) using the PuLP library.

    This class is responsible for building the linear programming model,
    including variables (routes), covering and fleet constraints,
    and solving the LP relaxation at each iteration of the column generation.
    """

    def __init__(self):
        self.mdl = pulp.LpProblem(sense=pulp.LpMinimize)
        self.mdl.setObjective(pulp.LpConstraintVar())
        self.solver = pulp.PULP_CBC_CMD(mip=False, msg=False)

        self.route_to_var: Dict[Route, pulp.LpVariable] = {}
        self.route_to_value: Dict[Route, float] = {}

        self.vertex_cover_constrs: Dict[vertex_id, pulp.LpConstraint] = {}
        self.vehicle_limit_constr: Optional[pulp.LpConstraint] = None

        self.obj: Optional[float] = None

    def init_cover_constr(self, u: vertex_id):
        """
        Initializes covering constraint for a customer.

        This constraint ensures that each customer is visited at least once.

        Args:
            u (vertex_id): The ID of the customer vertex for which to create the constraint.
        """

        constr = pulp.LpConstraint(-1, sense=pulp.LpConstraintGE, rhs=0, name=f"cover({u})")
        self.vertex_cover_constrs[u] = constr
        self.mdl.addConstraint(constr)

    def init_vehicle_limit_constr(self, vehicle_limit: int):
        """
        Initializes the constraint that limits the total number of vehicles used.

        Args:
            vehicle_limit (int): The maximum number of available vehicles.
        """

        constr = pulp.LpConstraint(-vehicle_limit, sense=pulp.LpConstraintLE, rhs=0, name="vehicle_limit")
        self.vehicle_limit_constr = constr
        self.mdl.addConstraint(constr)

    def add_route(self, route: Route):
        """
        Adds a new column (route) to the Restricted Master Problem.

        This involves creating a new PuLP variable for the route and updating the
        objective function and all relevant constraints (covering and vehicle limit)
        with this new variable.

        Args:
            route (Route): The route object to be added as a new column.
        """

        var = pulp.LpVariable(
            name=str(len(self.route_to_var)),
            lowBound=0,
            upBound=1,
            cat=pulp.LpContinuous
        )
        self.route_to_var[route] = var

        # Update objective and constraints
        self.mdl.objective.addInPlace(route.cost * var)
        for u in route.path:
            if u not in self.vertex_cover_constrs:
                continue
            constr = self.vertex_cover_constrs[u]
            constr.addInPlace(route.visit_count.get(u, 0) * var)
        if route.source != RouteSource.ARTIFICIAL:
            self.vehicle_limit_constr.addInPlace(var)

    def solve(self):
        """
        Solves the current LP relaxation of the RMP.

        After solving, it updates the model's objective value and the values of all
        route variables (columns).

        Raises:
            RuntimeError: If the LP model is not solved to optimality.
        """

        status = self.mdl.solve(self.solver)
        if status != pulp.LpStatusOptimal:
            raise RuntimeError("LP model is not solved to optimality")
        
        # Update objective and value of variables
        self.obj = self.mdl.objective.value()
        for route, var in self.route_to_var.items():
            self.route_to_value[route] = var.value()

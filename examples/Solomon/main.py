from cgvrp.vrptw.problem import Problem, VertexData, EdgeData
from cgvrp.vrptw.cg import ColumnGeneration

import sys
from pathlib import Path
import json

def load_solomon_instance(filepath):
    with open(filepath, "r") as f:
        data = json.load(f)
    problem = Problem()
    problem.vehicle_num = data["vehicle_number"]
    problem.vehicle_capacity = data["vehicle_capacity"]

    # Add vertices (customers)
    for customer in data["customers"]:
        vertex_data = VertexData(
            x=customer["x"],
            y=customer["y"],
            demand=customer["demand"],
            ready_time=customer["ready_time"],
            due_time=customer["due_time"],
            service_time=customer["service_time"]
        )
        is_depot = False
        if customer["id"] == 0:
            is_depot = True
        problem.add_vertex(customer["id"], vertex_data, is_depot)

    # Add edges
    for u, u_data in problem.vertices():
        for v, v_data in problem.vertices():
            if u == v:
                continue
            distance = ((u_data["x"] - v_data["x"])**2 + (u_data["y"] - v_data["y"])**2)**0.5
            edge_data = EdgeData(cost=distance, time=distance) # travel time = distance in Solomon dataset
            problem.add_edge(u, v, edge_data)

    return problem 

def main():
    if len(sys.argv) != 2:
        print("Usage: python main.py <instance_name>",
              "Example: python main.py c101", sep="\n")
        sys.exit(1)

    instance_name = sys.argv[1]
    filepath = Path(__file__).parent / "data" / f"{instance_name}.json"
    problem = load_solomon_instance(filepath)

    cg = ColumnGeneration(problem)
    result = cg.solve(gap_limit=1e-2, pricing_method="pulsing")
    # result = cg.solve(gap_limit=1e-2, pricing_method="labeling") # Use alternative pricing method

    print("\n" + "=" * 50)
    print(f"cost: {result.cost:.2f}")
    print(f"bound: {result.bound:.2f}")
    print(f"gap: {result.gap*100:.2f}%")

    print("\nRoutes:")
    for i, (route, lp_value) in enumerate(result.routes):
        print(f"  Route {i+1}: {route.path} (lp value: {lp_value:.4f})")

if __name__ == "__main__":
    main()
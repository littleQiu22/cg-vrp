import networkx as nx

from dataclasses import dataclass
from typing import Union, Optional

vertex_id = Union[str, int]

@dataclass
class VertexData:
    """
    Represents the attributes of a vertex (customer or depot) in the VRP.

    Attributes:
        x (float): The x-coordinate of the vertex.
        y (float): The y-coordinate of the vertex.
        demand (float): The product demand at the vertex.
        ready_time (float): The beginning of the time window in which service can start.
        due_time (float): The end of the time window by which service must start.
        service_time (float): The time required to perform the service at the vertex.
    """

    x: float
    y: float
    demand: float
    ready_time: float
    due_time: float
    service_time: float

@dataclass
class EdgeData:
    """
    Represents the attributes of an edge in the VRP graph.
    Attributes:
        cost (float): The cost of traversing the edge.
        time (float): The time required to traverse the edge.
    """

    cost: float
    time: float

class Problem:
    """
    Defines an instance of the Vehicle Routing Problem with Time Window (VRPTW).

    This class holds the problem's graph using NetworkX, as well as global
    information like the number of vehicles, their capacity, and the depot id.
    """

    def __init__(self):
        self.graph = nx.DiGraph()

        self.vehicle_num = 0
        self.vehicle_capacity = 0
        self.depot: Optional[vertex_id] = None

    def add_vertex(self, id: vertex_id, data: VertexData, is_depot: bool = False):
        """
        Adds a vertex to the problem graph.

        Args:
            id (vertex_id): The unique identifier for the vertex.
            data (VertexData): An object containing the vertex's attributes.
            is_depot (bool): Set to True if this vertex is the depot.
        """

        self.graph.add_node(
            id,
            x=data.x,
            y=data.y,
            demand=data.demand,
            ready_time=data.ready_time,
            due_time=data.due_time,
            service_time=data.service_time
        )
        if is_depot:
            self.depot = id

    def add_edge(self, u: vertex_id, v: vertex_id, data: EdgeData):
        """
        Adds a directed edge between two vertices in the graph.

        Args:
            u (vertex_id): The source vertex.
            v (vertex_id): The destination vertex.
            data (EdgeData): An object containing the edge's attributes.
        """

        self.graph.add_edge(
            u, v,
            cost=data.cost,
            time=data.time 
        )

    def vertices(self, vertex: vertex_id=None):
        """
        Gets vertex data from the graph.

        If no vertex is specified, acts as a generator yielding all vertices and their data.
        If a vertex is specified, returns the data for that particular vertex.

        Args:
            vertex (vertex_id, optional): The ID of the vertex to retrieve. Defaults to None.

        Yields:
            Tuple[vertex_id, dict]: The vertex ID and its attribute dictionary if no vertex is specified.

        Returns:
            dict: The attribute dictionary of the specified vertex.
        """
        if vertex is None:
            for u, data in self.graph.nodes(data=True):
                yield u, data
        else:
            return self.graph.nodes[vertex]

    def edges(self, prev: vertex_id=None, succ: vertex_id=None):
        """
        Gets edge data from the graph.

        If no edges are specified, acts as a generator yielding all edges and their data.
        If source and successor vertices are specified, returns the data for that particular edge.

        Args:
            prev (vertex_id, optional): The source vertex. Defaults to None.
            succ (vertex_id, optional): The destination vertex. Defaults to None.

        Yields:
            Tuple[vertex_id, vertex_id, dict]: The source, destination, and data of the edge.

        Returns:
            dict: The attribute dictionary of the specified edge.
        """
        if prev is None or succ is None:
            for u, v, edge in self.graph.edges(data=True):
                yield u, v, edge
        else:
            return self.graph.edges[prev, succ]
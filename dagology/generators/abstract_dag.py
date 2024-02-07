import networkx as nx
import numpy as np

__author__ = "\n".join(["Nicolas Kozak (nicolas.kozak5@gmail.com)"])
__all__ = ['AbstractDAG']

import dagology as dag

class AbstractDAG:
    def __init__(self):
        self.G = nx.DiGraph()
        self.N = 0
        self.topological_order = None

    def _add_nodes(self, coordinates):
        for i, pos in enumerate(coordinates):
            self.G.add_node(self.N + i, position=tuple(pos))
        self.N += len(coordinates)

    def _get_edges_from_neighbours(self, i, neighbours, weighted=False, weights=None):
        if weighted:
            return ((i, j, w) for j, w in zip(neighbours, weights))
        else:
            return ((i, j) for j in neighbours)
    
    def _add_edges_to_graph(self, edges, weighted=False):
        if weighted:
            self.G.add_weighted_edges_from(edges)
        else:
            self.G.add_edges_from(edges)

    def generate_graph(self, R, p=1.0, sorted=True, weighted=False):
        self._add_nodes(R)

        edges = []
        for i in range(self.N - 1):
            j_start = i + 1 if sorted else 0
            i_weights = self._calculate_weights(R[i], R[j_start:])
            i_nbrs = self._is_neighbor(R[i], R[j_start:], i_weights)
            neighbours = np.where(i_nbrs)[0] + j_start
            edges.extend(self._get_edges_from_neighbours(i, neighbours, weighted, i_weights[i_nbrs]))
            
        if p < 1:
            mask = np.random.rand(len(edges)) > p
            edges = np.array(edges)[mask].tolist()

        self._add_edges_to_graph(edges, weighted)
        if sorted:
            self.topological_order = range(self.N)
        return self.G
    
    def traverse_path(self, update_path, traversal_type='forward'):
        """
        Traverse the DAG either forward or backward in topological order.
        Returns a path according to an update rule.
        
        :param update_path: The function defining the update rule.
        :param traversal_type: Type of traversal: 'forward' or 'backward'. Defaults to 'forward'.

        Example update_path rule for computing shortest path:
        >>>
        def update_path(u, weights, d, p):
            for v, w in weights.items():
                new_d = d[u] + w
                if d[v] > new_d:
                    d[v], p[v] = new_d, u
        >>>
        """
        if traversal_type not in ['forward', 'backward']:
            raise ValueError("Invalid traversal type. Type should be 'forward' or 'backward'.")

        if not self.topological_order:
            raise ValueError("Provide a valid topological order")

        distances = [float("inf")] * self.N
        predecessors = {}

        source, target = self.topological_order[0], self.topological_order[-1]

        if traversal_type == 'forward':
            distances[source] = 0
            traversal_order = self.topological_order
            get_neighbors = self.G.neighbors
        else:
            distances[target] = 0
            traversal_order = reversed(self.topological_order)
            get_neighbors = self.G.predecessors

        edge_weights = nx.get_edge_attributes(self.G, 'weight')

        for node in traversal_order:
            neighbors = list(get_neighbors(node))
            if traversal_type == 'forward':
                weights = {v: edge_weights[(node, v)] for v in neighbors}
            else:
                weights = {v: edge_weights[(v, node)] for v in neighbors}
            update_path(node, weights, distances, predecessors)

        path = []
        current_node = target if traversal_type == 'forward' else source
        while current_node is not None:
            path.append(current_node)
            current_node = predecessors.get(current_node)
        path.reverse()

        return path

    def _is_neighbor(self, x, y, weights):
        raise NotImplementedError("Subclasses must implement this method")

    def _calculate_weights(self, x, y):
        raise NotImplementedError("Subclasses must implement this method")

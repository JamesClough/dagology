import networkx as nx
import numpy as np

__author__ = "\n".join(["Nicolas Kozak (nicolas.kozak5@gmail.com)"])
__all__ = ['AbstractDAG']

import dagology as dag

class AbstractDAG:
    def __init__(self):
        self.G = nx.DiGraph()
        self.N = 0

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
        return self.G

    def _is_neighbor(self, x, y, weights):
        raise NotImplementedError("Subclasses must implement this method")

    def _calculate_weights(self, x, y):
        raise NotImplementedError("Subclasses must implement this method")

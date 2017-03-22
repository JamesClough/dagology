from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true
from nose.tools import assert_almost_equals

import networkx as nx
import numpy as np
import dagology as dag


class TestEmbed(object):
    """ Unit tests for embedding"""

    def test_shape(self):
        G = nx.DiGraph()
        G.add_edges_from([[0,1], [0,2], [0,3],
                          [1,4], [2,4],
                          [1,5], [3,5],
                          [2,6], [3,6],
                          [4,7], [5,7], [6,7]])
        for D in [1,2,3]:
            X, nodes = dag.minkowski_embed(G, D)
            assert_equal(X.shape, (8,D))

    def test_line(self):
        G = nx.path_graph(5, create_using=nx.DiGraph())
        X, nodes = dag.minkowski_embed(G, 2)
        if X[0,0] > X[0,-1]:
            X *= -1
        correct_coords = ((-2, 0),
                          (-1, 0),
                          (0, 0),
                          (1, 0),
                          (2, 0))
        for i in range(5):
            for j in [0,1]:
                assert_almost_equals(X[i,j], correct_coords[i][j])
                
    def test_dimension(self):
        """ Check first embedding coordinates not changed by adding more
        dimensions"""
        G = nx.DiGraph()
        for edge in nx.petersen_graph().edges():
            if edge[0] < edge[1]:
                G.add_edge(edge[0], edge[1])
        X_2, nodes = dag.minkowski_embed(G, 2)
        X_5, nodes = dag.minkowski_embed(G, 5)
        for i in range(10):
            for j in range(2):
                assert_equal(X_2[i,j], X_5[i,j])
        

from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true

import networkx as nx
import numpy as np
import dagology as dag


class TestRandomDag(object):
    """ Unit tests for the Newman Karrer random DAG model"""
    def test_example(self):
        degree_sequence = [[0, 2], [1, 1], [2,0]]
        G = dag.random_dag(degree_sequence)
        assert_equal(G.number_of_nodes(), 3)
        assert_equal(G.number_of_edges(), 3)
        
    def test_multiedge(self):
        degree_sequence = [[0, 2], [0, 2], [2, 0], [2, 0]]
        G = dag.random_dag(degree_sequence)
        assert_equal(G.number_of_nodes(), 4)
        assert_equal(G.number_of_edges(), 4)
        assert_true(G.has_edge(0, 2))
        assert_true(G.has_edge(0, 3))        
        

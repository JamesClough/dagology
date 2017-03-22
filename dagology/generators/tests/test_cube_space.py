from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true

import networkx as nx
import dagology as dag

class TestCubeSpace(object):
    """ Unit tests for the cube space model"""
    def test_number_of_nodes(self):
        N = 80
        for p in [0.0, 0.1, 0.9, 1.0]:
            G = dag.cube_space_graph(N, 5, p)
            assert_equal(G.number_of_nodes(), N)

    def test_number_of_edges(self):
        N = 100
        G = dag.cube_space_graph(N, 2, 0.)
        assert_true(G.number_of_edges() == 0)

        G = dag.cube_space_graph(N, 1, 1.)
        assert_equal(G.number_of_edges(), (N*(N-1)/2))

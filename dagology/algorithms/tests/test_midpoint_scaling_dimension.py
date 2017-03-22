from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true

import networkx as nx
import numpy as np
import dagology as dag


class TestMPSD(object):
    """ Unit tests for midpoint scaling dimension"""

    def test_line(self):
        G = nx.DiGraph()
        G.add_edges_from([[x, x + 1] for x in range(10)])
        assert_equal(dag.mpsd(G), 1.0)

    def test_diamond(self):
        G = nx.DiGraph()
        G.add_edges_from([[0, 1], [0, 2], [1, 3], [2, 3],
                          [3, 4], [3, 5], [3, 6],
                          [4, 7], [5, 7], [6, 7],
                          [1, 8], [8, 9], [9, 10], [10, 11], [11, 4],
                          [2, 12], [12, 13], [13, 14], [14, 15], [15, 6]])
        assert_equal(dag.mpsd(G), 2.0)

    def test_empty(self):
        G = nx.DiGraph()
        assert_equal(dag.mpsd(G), 0)

    def test_unconnected(self):
        G = nx.DiGraph()
        G.add_nodes_from([1, 2, 3])
        assert_equal(dag.mpsd(G), 0)

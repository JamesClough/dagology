from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true
from nose.tools import assert_almost_equals

import networkx as nx
import numpy as np
import dagology as dag


class TestMMD(object):
    """ Unit tests for Myrheim-Meyer dimension"""

    def test_line(self):
        G = nx.DiGraph()
        G.add_edges_from([[x, x + 1] for x in range(10)])
        assert_almost_equals(dag.mmd(G), 1.16, 5)  # small N correction

    def test_empty(self):
        G = nx.DiGraph()
        assert_equal(dag.mmd(G), 0)

    def test_unconnected(self):
        G = nx.DiGraph()
        G.add_nodes_from([1, 2, 3])
        assert_equal(dag.mmd(G), 0)


class TestMMDFormula(object):
    """ Unit tests for the MMD formula"""

    def test_bounds(self):
        assert_raises(ValueError, dag.mmd_formula, 0.5, 2)
        assert_raises(ValueError, dag.mmd_formula, -1, 3)

    def test_k2(self):
        assert_equal(dag.mmd_formula(2, 2), 0.25)
        assert_almost_equals(dag.mmd_formula(3, 2), (4. / 35.), 5)
        assert_almost_equals(dag.mmd_formula(4, 2), (1. / 20.), 5)

    def test_k3(self):
        assert_almost_equals(dag.mmd_formula(2, 3), (1. / 36.), 5)
        assert_almost_equals(dag.mmd_formula(3, 3), (2. / 525.), 5)
        assert_almost_equals(dag.mmd_formula(4, 3), (1. / 2100.), 6)

    def test_k4(self):
        assert_almost_equals(dag.mmd_formula(2, 4), (1. / 576.), 6)
        assert_almost_equals(dag.mmd_formula(3, 4), (4. / 75075.), 7)
        assert_almost_equals(dag.mmd_formula(4, 4), (1. / 705600.), 8)

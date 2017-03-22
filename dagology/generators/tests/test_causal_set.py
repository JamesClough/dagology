from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true

import networkx as nx
import numpy as np
import dagology as dag


class TestCausalSet(object):
    """ Unit tests for the causal set model"""
    def test_number_of_nodes(self):
        N = 80
        D = 2
        R = np.random.random((N, D)) # probably shouldn't test stochastically...
        for p in [0.0, 0.1, 0.9, 1.0]:
            G = dag.causal_set_graph(R, p)
            assert_equal(G.number_of_nodes(), N)

    def test_number_of_edges(self):
        N = 80
        D = 2
        R = np.random.random((N, D))
        G = dag.causal_set_graph(R, 0.)
        assert_true(G.number_of_edges() == 0)

    def periodicity(self):
        R = np.array([[0., 0.],
                      [1., -0.8],
                      [1., 0.8],
                      [2., 1.6]])
        G_boundary = dag.causal_set_graph(R)
        G_periodic = dag.causal_set_graph(R, period=2.)
        assert_equal(G_boundary.number_of_edges(), 4)
        assert_equal(G_periodic.number_of_edgeS(), 5)

class TestMinkowskiInterval(object):
    """ Unit tests for minkowski_interval"""
    def test_shape(self):
        N = 120
        D = 3
        R = dag.minkowski_interval(N, D)
        assert_equal(R.shape, (N,D))

    def test_fix_ends_true(self):
        N = 100
        D = 2
        R = dag.minkowski_interval(N, D, fix_ends=True)
        assert_equal(R[0, 0], 0.)
        assert_equal(R[0, 1], 0.5)
        assert_equal(R[1, 0], 1.)
        assert_equal(R[1, 1], 0.5)

    def test_fix_ends_false(self):
        N = 100
        D = 2
        R = dag.minkowski_interval(N, D, fix_ends=False)
        assert_true(0. < R[0, 0] < 1.)
        assert_true(R[0, 1] != 0.5)
        assert_true(0. < R[0, 0] < 1.)
        assert_true(R[0, 1] != 0.5)

class TestDeSitterInterval(object):
    """ Unit tests for de_sitter_interval"""
    def test_shape(self):
        N = 120
        D = 3
        R = dag.de_sitter_interval(N, D, 0.1)
        assert_equal(R.shape, (N,D))

    def test_fix_ends_false(self):
        N = 100
        D = 2
        R = dag.de_sitter_interval(N, D, 0.5, fix_ends=False)
        assert_true(0. < R[0, 0] < 1.)
        assert_true(R[0, 1] != 0.5)
        assert_true(0. < R[0, 0] < 1.)
        assert_true(R[0, 1] != 0.5)

    def test_fix_ends_true(self):
        N = 100
        D = 2
        R = dag.de_sitter_interval(N, D, 0.5, fix_ends=True)
        assert_equal(R[0, 0], 0.)
        assert_equal(R[0, 1], 0.)
        assert_equal(R[1, 0], 1.)
        assert_equal(R[1, 1], 0.)

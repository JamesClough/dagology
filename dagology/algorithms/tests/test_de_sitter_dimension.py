from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true
from nose.tools import assert_almost_equals

import networkx as nx
import numpy as np
import dagology as dag
from dagology.dagology.algorithms.de_sitter_dimension import C_1, C_2, C_3

class TestDeSitterFormula(object):
    """ Unit tests for de Sitter dimension
    
    Note - this should currently be considered as a work in progress
    For some areas of the T,d,K parameter space it seems the numerical
    optimisation struggles with some non-unique, or very close sets of
    chain values."""

    def test_equation_solver_1(self):
        allowed_error = 0.2
        T, d, K = 20., 2., 0.002
        chains = (C_1(T,d,K), C_2(T,d,K), C_3(T,d,K))
        T_, d_, K_, = dag.de_sitter_param_estimate(chains)
        assert_true(abs(T - T_)/T < allowed_error)
        assert_true(abs(d - d_)/d < allowed_error)
        assert_true(abs((K - K_)/K) < allowed_error)




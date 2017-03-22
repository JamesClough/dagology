from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true

import networkx as nx
import numpy as np
import dagology as dag

class TestLongestPathMatrix(object):
    """ Unit tests for longest_path_matrix function"""

    def test_line(self):
        G = nx.path_graph(10, create_using=nx.DiGraph())
        A = nx.adjacency_matrix(G)
        A = A.toarray()
        LP = dag.longest_path_matrix(A)
        for i in range(10):
            for j in range(10):
                assert_equal(LP[i,j], max(j-i, 0))

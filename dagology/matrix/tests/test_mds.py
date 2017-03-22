from nose.tools import assert_equal
from nose.tools import assert_almost_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true

import networkx as nx
import numpy as np
import dagology as dag


class TestMDS(object):

    def test_line(self):
        N = 10
        ds2 = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                ds2[i,j] = (i-j)**2
        X = dag.mds(ds2, 1)
        for i in range(N):    
            assert_almost_equal(X[i, 0], (N-1)/2. - i)
        

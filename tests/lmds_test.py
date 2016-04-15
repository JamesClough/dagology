""" Test Lorentzian Multidimensional Scaling"""

import dagology as dag
import numpy as np
import time

def test_minkowski_scatter():
    N = 1000
    D = 2
    R = dag.minkowski_interval(N, D)
    ds2 = dag.sq_separations(R, dag.minkowski)
    X = dag.LMDS(ds2, D, k=20, method='lorentzian')
    assert X.shape == (N, D)
    
    
def main():
    test_minkowski_scatter()
    
if __name__ == "__main__":
    main()

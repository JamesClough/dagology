""" Test Dimension measures"""
import dagology as dag
import numpy as np
import time

def test_causal_set():
    N, D = 500, 2
    t = time.time()
    A = dag.interval_causet(N, D, fix_ends=True)
    print 'Creating Causet with N=%s -- D=%s - %s sec' % (N, D, np.round(time.time() - t, 2))
    t = time.time()
    D_MMD = dag.myrheim_meyer(A, k=2)
    print 'Myrheim-Meyer - %s sec' % (np.round(time.time() - t, 2))
    t = time.time()
    D_MPSD = dag.midpoint_scaling(A)
    print 'Midpoint-Scaling - %s sec' % (np.round(time.time() - t, 2))
    print D_MMD, D_MPSD
    
def main():
    test_causal_set()
    
if __name__ == "__main__":
    main()

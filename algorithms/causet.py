""" Calculate Minkowski separations from adjacency matrix of a causal set or DAG"""
import numpy as np
import metrics
   
def causet_adj_matrix(S, R):
    """ Return causal set adjacency matrix A
    
        S: separations
        R: original coordinates"""
    N = S.shape[0]
    A = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            # check time ordering - A[i,j] is 1 if i is in the future of j
            if R[i,0] > R[j,0]:
                if S[i,j] < 0:
                    A[i,j] = 1.
    return A  

def interval_causet(N, D, fix_ends=True):
    """ Return N point, D dimensional Minkowski scattered causet adj-matrix"""
    import scatter
    import dagsep
    R = scatter.minkowski_interval(N, D, fix_ends)
    S = metrics.sq_separations(R, metrics.minkowski)
    A = causet_adj_matrix(S, R)
    return A            

if __name__ == "__main__":
    print __doc__
                

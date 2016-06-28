""" Create adjacency matrices for causal sets"""
import numpy as np
import metrics
import dagsep   
   
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

def minkowski_causet(N, D, fix_ends=True):
    """ Return N point, D dimensional Minkowski scattered causet adj-matrix
        A - adjacency matrix     NxN
        R - original coordinates NxD"""
    import scatter
    import dagsep
    R = scatter.minkowski_interval(N, D, fix_ends)
    S = metrics.sq_separations(R, metrics.minkowski)
    A = causet_adj_matrix(S, R)
    return A, R     
    
def de_sitter_causet(N, D, eta_0, eta_1):
    """ Return N point, D dimensional de Sitter scattered causet adj-matrix
        A - adjacency matrix     NxN
        R - original coordinates NxD"""
    import scatter
    import dagsep
    R = scatter.de_sitter_interval(N, D, eta_0, eta_1)
    S = metrics.sq_separations(R, metrics.de_sitter)
    A = causet_adj_matrix(S, R)
    return A, R       

def random_dag(N, k, tc=True):
    p = (1. * k)/N
    A = np.random.random((N,N))
    A[A>(1-p)] = 1
    A[A<(1-p)] = 0
    for i in range(N):
        A[i,i:] = 0
    if tc:
        A = dagsep.transitive_completion(A)
    return A

if __name__ == "__main__":
    print __doc__
                

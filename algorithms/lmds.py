""" Lorentzian Multidimensional scaling"""
import numpy as np
from scipy.spatial.distance import pdist, squareform
from numpy.linalg import eigh

################################################################################
# Classic MDS algorithm ########################################################
def J_matrix(N):
    """ Return NxN double-centering matrix"""
    J = np.identity(N) - (1./N)*np.ones((N, N))
    return J
    
def calc_double_centre(S):
    """ Return double-centered square distance matrix
    
    S -- squared distance matrix"""
    N, _ = S.shape
    assert N == _ , 'Distance matrix must be square'
    J = J_matrix(N)
    A = -0.5 * np.dot(np.dot(J, S), J)
    return A

def e_eval(E, U, D):
    """ Return largest D eigenvectors
    
    Euclidean eigenvalues - extracts D largest eigenvalues and eigenvectors
    E -- N array of eigenvalues
    U -- NxN array of eigenvectors
    D -- Dimension (integer)
    
    Returns:
    E_max -- Diagonal DxD array of sqrt of D largest eigenvalues
    U_max -- NxD array of corresponding eigenvectors"""
    E_max = np.diag(E[-1*D:])
    U_max = U[:, -1*D:]
    return E_max, U_max
    
def m_eval(E, U, D):
    """ Return largest D Lorentzian eigenvalues
    
    Minkowski eigenvalues - extracts D-1 largest eigenvalues and eigenvectors
                            and the largest negative eigenvalue and vector
    E -- N array of eigenvalues
    U -- NxN array of eigenvectors

    Returns:
    E_max -- Diagonal DxD array of sqrt of eigenvalues
    U_max -- NxD array of corresponding eigenvectors"""
    E_max = np.diag(np.concatenate((E[0:1]*-1., E[-1*(D-1):])))
    U_max = np.concatenate((U[:, 0:1], U[:, -1*(D-1):]), axis=1)
    return E_max, U_max

def MDS(ds2, D, method='euclidean'):
    """ Classic MDS algorithm
    
    Allowed methods - 
    -- euclidean - classic MDS
    -- lorentzian - the largest negative eigenvalue is used, and the D-1 largest positive"""
    N, _ = ds2.shape
    assert N == _
    A = calc_double_centre(ds2)
    E_, U_ = eigh(A)
    if method == 'euclidean':
        E, U = e_eval(E_, U_, D)
    elif method == 'lorentzian':
        E, U = m_eval(E_, U_, D)
    X = np.dot(U, np.sqrt(E))
    return X
       
################################################################################
# Landmark MDS
def select_landmarks(ds2, k):
    """ Choose k landmark points 
    
    Select k landmark points and shuffle ds2 matrix so that these k points
    are the first k in the index"""
    # do random for now
    return ds2 

def LMDS(ds2, D, k=None, method='euclidean'):
    """ Landmark MDS algorithm 
    
    Instead of fitting every point to every other point, first fit k landmark
    points to each other, and then fit the remaining N-k to those landmarks
    If the landmarks are well distributed then this is a good approximation
    of MDS on all N points.
    JC - should be O(k*n) in speed """
    N, _ = ds2.shape
    assert N == _
    if not k:
        k = min(10, D*2)
    ds2 = select_landmarks(ds2, k)
    ds2_L = ds2[:k, :k]
    ds2_rest = ds2[k:, :k]
    
    A_L = calc_double_centre(ds2_L)
    E_, U_ = eigh(A_L)
    if method == 'euclidean':
        E, U = e_eval(E_, U_, D)
    elif method == 'lorentzian':
        E, U = m_eval(E_, U_, D)
    
    X_L = np.dot(U, np.sqrt(E))
    mean_rest = np.mean(ds2_rest, axis=0)
    L = U * 1./np.sqrt(np.diag(E))
    X = np.dot((mean_rest - ds2_rest), L) * 0.5
    # temporary hack I'm sorry
    if method=='lorentzian':
        X[:, 0] *= -1
    return np.concatenate((X_L, X), axis=0)

################################################################################
# Pivot MDS 
def pivot_double_centre(C):
    """ C - Nxk matrix of square distances from pivots to points
        Returns double-centered distance matrix of same"""
    N, k = S.shape
    assert N > k
    A = np.zeros((N, k))
    for i in range(N):
        for j in range(k):
            if i==j:
                A[i,j] += 1
            A[i,j] -= np.sum(C[:,j]) * (1./N)
            A[i,j] -= np.sum(C[i,:]) * (1./k)
            A[i,j] += np.sum(C) * (1./(k*N))
    return A*-0.5

def PMDS(ds2, D, k=None, method='euclidean'):
    """ Pivot MDS algorithm
    
    NOT FULLY IMPLEMENTED AND TESTED YET
    """
    N, _ = ds2.shape
    assert N == _
    if not k:
        k = min(10, D*2)
    # going to guess at how this works first - JC
    ds2 = select_landmarks(ds2, k)
    ds2_P = ds2[:, :k] # this is the matrix C in Brandes paper
    ds2_rest = ds2[k:, k:]
    A_L = pivot_double_centre(ds2_P)
    CC = np.dot(A_L, A_L.trans)
    E_, U_ = eigh(CC)
    if method == 'euclidean':
        E, U = e_eval(E_, U_, D)
    elif method == 'lorentzian':
        E, U = m_eval(E_, U_, D)
    X_P = np.dot(U, np.sqrt(E))
    mean_rest = np.mean(ds2_rest, axis=0)
    P = U * 1./np.sqrt(np.diag(E))
    X = np.dot((mean_rest - ds2_rest), P) * 0.5
    # temporary hack I'm sorry
    if method=='lorentzian':
        X[:, 0] *= -1
    return np.concatenate((X_P, X), axis=0)

if __name__ == "__main__":
    print __doc__                



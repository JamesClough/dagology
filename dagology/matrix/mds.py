""" Lorentzian Multidimensional scaling"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)"])

import numpy as np
from numpy.linalg import eigh

__all__ = ['mds']


def J_matrix(N):
    """ Return NxN double-centering matrix"""
    J = np.identity(N) - (1./N)*np.ones((N, N))
    return J
    
def calc_double_centre(S):
    """ Return double-centered square distance matrix
    
    Parameters
    ----------
    
    S : numpy array
        matrix of square distances
        
    Returns
    -------
    
    Double centered square distance matrix"""
    
    N, _ = S.shape
    assert N == _ , 'Distance matrix must be square'
    J = J_matrix(N)
    A = -0.5 * np.dot(np.dot(J, S), J)
    return A

def e_eval(E, U, D):
    """ Return largest D eigvenvalues and corresponding eigenvectors
    
    Euclidean eigenvalues - extracts D largest eigenvalues and eigenvectors
    E -- N array of eigenvalues
    U -- NxN array of eigenvectors
    D -- Dimension (integer)
    
    Returns
    -------
    
    E_max : Diagonal DxD array of sqrt of D largest eigenvalues
    U_max : NxD array of corresponding eigenvectors"""
    E_max = np.diag(E[-1*D:])
    U_max = U[:, -1*D:]
    return E_max, U_max
    
def m_eval(E, U, D):
    """ Return largest D Lorentzian eigvenvalues and corresponding eigenvectors
    
    Minkowski eigenvalues - extracts D-1 largest eigenvalues and eigenvectors
                            and the largest negative eigenvalue and vector
    E -- N array of eigenvalues
    U -- NxN array of eigenvectors
    
    Returns
    -------
    
    E_max : Diagonal DxD array of sqrt of eigenvalues
    U_max : NxD array of corresponding eigenvectors"""
    
    E_max = np.diag(np.concatenate((E[:1]*-1., E[:-1*D:-1])))
    U_max = np.concatenate((U[:, :1], U[:, :-1*D:-1]), axis=1)
    return E_max, U_max

def mds(ds2, D, method='euclidean'):
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
    else:
        assert False, 'method must be either euclidean or lorentzian'
    X = np.dot(U, np.sqrt(E))
    return X
       


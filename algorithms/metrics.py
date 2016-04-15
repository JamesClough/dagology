""" Provide various metrics for calculating distances between two points 
General form - metric(x, y) returns real number distance between points x and y
"""
import numpy as np

"""
################################################################################
Riemannian
################################################################################
"""
def euclidean(x, y):
    """Calculcate Euclidean distance between x and y """
    assert len(x)==len(y), 'ERROR - vectors in Euclidean metric have different lengths'
    return np.sum((np.array(x) - np.array(y))*(np.array(x) - np.array(y)))

"""
################################################################################
Lorentzian
################################################################################
"""
def minkowski(x, y):
    """Calculate Minkowski separation between x and y using -++...+ convention"""
    assert len(x)==len(y), 'ERROR - vectors in Minkowski metric have different lengths'
    dt = x[0] - y[0]
    dt2 = dt * dt
    dx = np.array([x[i] - y[i] for i in range(1, len(x))])
    dx2 = dx * dx
    dx2sum = sum(dx2)
    return dx2sum - dt2 
      
def minkowski_periodic(x, y, L=[]):
    """Calculate Minkowski separation between x and y using -++...+ convention
       Periodic boundary conditions in spatial coordinates are given in L
       If len(L) < D-1 then assume no boundary on other spatial dimensions"""
    D = len(x)
    assert len(y) == D, 'ERROR - vectors in Minkowski metric have different lengths'
    # fill in remaining box dimensions with None
    # JC - we can probably make this more efficient later by vectorising the
    #      operation if the if L_d line
    while len(L) < (D-1):
        L.append(None)
    dt = x[0] - y[0]
    dt2 = dt * dt
    ds2 = -1 * dt2
    for d in range(1, D):
        L_d = L[d-1] # the first dimension in the coords is time so exclude it
        if L_d:
            dx2 = min((x[d]-y[d])**2, (x[d]-y[d]+L_d)**2, (x[d]-y[d]-L_d)**2)
        else:
            dx2 = (x[d]-y[d])**2
        dz2 += dx2
    return ds2

def sq_separations(R, metric):
    """ Return square array of squared Minkowski separations from coords R
        c - speed of light"""
    S = np.zeros((R.shape[0], R.shape[0]))
    for i in range(R.shape[0]):
        for j in range(R.shape[0]):
            S[i, j] = metric(R[i,:], R[j,:])
    return S
    
if __name__ == "__main__":
    print __doc__
        
        
        
        
        
    

""" Generate uniform random coords in interval of space or spacetime"""
import numpy as np
import metrics

################################################################################
# MINKOWSKI
def minkowski_interval_scatter(N, D, fix_ends=True):
    """ Scatter N points in a D dimensional interval in Minkowski space
    
    Throw points into a unit box rejecting those outside the interval
    Repeat until N points have been reached
    Note that this is inefficient for large D"""
    R = np.random.random((N, D))
    if fix_ends:
        R[0, 0] = 0.
        R[1, 0] = 1.
        R[0:2, 1:] = 0.5
    a = R[0, :]
    b = R[1, :]
    for i in range(2, N):
        while (metrics.minkowski(a, R[i, :]) > 0) or ((metrics.minkowski(R[i, :], b) > 0)):
            R[i, :] = np.random.random(D)
    return R

def minkowski_interval_map(N, D, fix_ends=True):
    """ Scatter N points in a D dimensional interval in Minkowski space
    
    Build Minkowski interval in `clever' way by mapping [0,1]^D to 
    the correct spacetime coords
    """
    assert False, 'ERROR - minkowski_interval_map not implemented yet'
    
def minkowski_interval(N, D, fix_ends=True, method='scatter'):
    """ Scatter N points in a D dimensional interval in Minkowski space
    
    Available methods are:
    scatter -- place points in unit cube and check they lie within
               the appropriate interval, keeping those that do.
               This is slow for large D - something like 2^D slowdown
               
    map -- map D unit cube to the relevant interval respecting volume elements
           not yet implemented
    """
    if method == 'scatter':
        return minkowski_interval_scatter(N, D, fix_ends=True)
    elif method == 'map':
        return minkowski_interval_map(N, D, fix_ends=True)
    else:
        assert False, 'Invalid method %s given to minkowski_interval' % method

################################################################################
# DE SITTER
# NOT IMPLEMENTED YET
def de_sitter_interval(N, D, a, fix_ends=True):
    assert False, 'ERROR - de sitter interval not implemented yet'
    
    
if __name__ == "__main__":
    print __doc__

         

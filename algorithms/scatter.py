""" Generate uniform random coords in interval of space or spacetime"""
import numpy as np
import metrics

################################################################################
# SPHERE
def sphere_surface_cartesian(N, D):
    """ Generate N points uniformly sampled from surface of a D-sphere
    
    Return Cartesian coordinates
    Using normal distributions as multivariate normal is spherically symmetric
    """
    R = np.random.randn(N, D+1)
    R_sq = R*R
    R_sq_sum = np.sqrt(np.sum(R_sq, axis=1))
    R_norm = R_sq_sum.reshape(N, 1)
    return R/R_norm

def sphere_surface_angular(N, D):
    """ Generate N points uniformly sampled from surface of a D-sphere"""
    X = sphere_surface_cartesian(N, D)
    R = np.zeros((N, D))
    for i in range(N):
        R[i,:] += metrics.cartesian_to_angular(X[i,:])
    return R

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
# HYPERBOLIC
def hyperbolic_disk(N, R, a=1.):
    """ Scatter N points in a 2 dimensional hyperbolic manifold with curvature a
    
    The points are scattered uniformly with inside a disk of radius R 
    We are using the native representation, where polar coordinate r
    is the hyperbolic distance to the origin"""
    X = np.random.rand(N, 2)
    X[:,1] *= (2.* np.pi)
    A_R = np.cosh(R * a) - 1.
    X[:,0] = np.arccosh((X[:,0] * A_R) + 1.) / a
    return X
    
################################################################################
# DE SITTER
def de_sitter_interval(N, D, eta_0, eta_1, method='scatter'):
    """ Scatter N points in a D dimensional interval in de Sitter space
    
    Available methods are:
    scatter -- place points in unit cube and check they lie within
               the appropriate interval, keeping those that do.
               This is slow for large D - something like 2^D slowdown
               
    map -- map D unit cube to the relevant interval respecting volume elements
           not yet implemented
    """
    if method == 'scatter':
        return de_sitter_interval_scatter(N, D, eta_0, eta_1)
    elif method == 'map':
        return de_sitter_interval_map(N, D, eta_0, eta_1)
    else:
        assert False, 'Invalid method %s given to de_sitter_interval' % method    
    
def de_sitter_interval_scatter(N, D, eta_0, eta_1):
    """ Scatter N points in a D dimensional de Sitter space
    
    We are using conformal coordinates: eta (conformal time)
                                        theta (spatial coords)
    
    which works out the same as flat spacetime except with the points uniformly
    scattered in tan(eta)
    
    arguments: N -- number of points
               D -- dimension of spacetime
           eta_0 -- starting conformal time - choose eta_0 > 0
           eta_1 -- ending conformal time   - choose eta_1 < pi/2 
           
    JC - this is a sprinkling method so can be inefficient if many points are
         thrown away - this occurs when D is large, or when an eta_1 is close
         to pi/2
         
    for now, assume we want to fix the endpoints
    """

    u_0, u_1 = np.tan(eta_0), np.tan(eta_1)
    d_u = u_1 - u_0
    R = np.zeros((N, D))
    R[0,0] += eta_0 
    R[1,0] += eta_1
    i = 2
    while i < N:
        u = (np.random.random() * d_u) + u_0
        eta_u = np.arctan(u)
        theta = sphere_surface_angular(1, D-1)[0]
        x = np.concatenate(([eta_u], theta))
        if metrics.de_sitter(x, R[0,:]) < 0:
            if metrics.de_sitter(x, R[1,:]) < 0:
                R[i,:] += x
                i += 1
    return np.array(R)

    
if __name__ == "__main__":
    print __doc__

         

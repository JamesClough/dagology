"""
Models for causal set graphs.
"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)", "Nicolas Kozak (nicolas.kozak5@gmail.com)"])

import networkx as nx
import numpy as np

import dagology as dag
from .abstract_dag import AbstractDAG

__all__ = ['CausalSetGraph',
           'minkowski_interval',
           'de_sitter_interval']

class CausalSetGraph(AbstractDAG):
    def _is_neighbor(self, x, y, weights):
        return np.logical_and(weights > 0, x[0] < y[:, 0])

    def _calculate_weights(self, x, y):
        return dag.minkowski(x, y)


def minkowski_interval_scatter(N, D, fix_ends=True):
    """ Scatter N points in a D dimensional interval in Minkowski space.
    Start point has coordinates: (0, 0.5, 0.5, ..., 0.5).
    End point has coordinates: (1, 0.5, 0.5, ..., 0.5).
    +--- metric is assumed.

    Parameters
    ----------

    N - number of points
    D - dimension of spacetime
    fix_ends - if True, have points at start and end of interval

    Notes
    -----

    Throw points into a unit box rejecting those outside the interval
    Repeat until N points have been reached
    Note that this is inefficient for large D"""
    
    R = np.random.random((N, D))
    a, b = 0.5*np.ones(D), 0.5*np.ones(D)
    a[0], b[0] = 0., 1.

    if fix_ends:
        R[0], R[N-1] = a, b
    for i in range(1, N-1) if fix_ends else range(N):
        while dag.minkowski(a, R[i]) < 0 or dag.minkowski(R[i], b) < 0:
            R[i] = np.random.random(D)
    return R


def minkowski_interval_map(N, D, fix_ends=True):
    """ Scatter N points in a D dimensional interval in Minkowski space

    Build Minkowski interval in `clever' way by mapping [0,1]^D to
    the correct spacetime coords
    """
    assert False, 'ERROR - minkowski_interval_map not implemented yet'


def minkowski_interval(N, D, fix_ends=True, method='scatter', sorted=True):
    """ Scatter N points in a D dimensional interval in Minkowski space

    Available methods are:
    scatter -- place points in unit cube and check they lie within
               the appropriate interval, keeping those that do.
               This is slow for large D - something like 2^D slowdown

    map -- map D unit cube to the relevant interval respecting volume elements
           not yet implemented
           
    Parameters
    ----------
    sorted - boolean - If true, sorts the array's zeroth component in ascending order
    """
    methods = {'scatter': minkowski_interval_scatter,
               'map': minkowski_interval_map}
    
    if method not in methods:
        raise ValueError(f"Invalid method {method} given to minkowski_interval")

    R = methods[method](N, D, fix_ends)
    if sorted:
        # Sort by increasing time component
        R = R[np.argsort(R[:, 0])]
    return R


def sphere_surface_cartesian(N, D):
    """ Generate N points uniformly sampled from surface of a D-sphere

    Return Cartesian coordinates
    Using normal distributions as multivariate normal is spherically symmetric
    """
    R = np.random.randn(N, D + 1)
    R_sq = R * R
    R_sq_sum = np.sqrt(np.sum(R_sq, axis=1))
    R_norm = R_sq_sum.reshape(N, 1)
    return R / R_norm


def sphere_surface_angular(N, D):
    """ Generate N points uniformly sampled from surface of a D-sphere"""
    X = sphere_surface_cartesian(N, D)
    R = np.zeros((N, D))
    for i in range(N):
        R[i, :] += dag.cartesian_to_angular(X[i, :])
    return R


def hyperbolic_disk(N, R, a=1.):
    """ Scatter N points in a 2 dimensional hyperbolic manifold with curvature a

    The points are scattered uniformly with inside a disk of radius R
    We are using the native representation, where polar coordinate r
    is the hyperbolic distance to the origin"""
    X = np.random.rand(N, 2)
    X[:, 1] *= (2. * np.pi)
    A_R = np.cosh(R * a) - 1.
    X[:, 0] = np.arccosh((X[:, 0] * A_R) + 1.) / a
    return X

def de_sitter_interval(N, D, KT2, fix_ends=False, method='scatter'):
    if method == 'scatter':
        return de_sitter_interval_scatter(N, D, KT2, fix_ends)
    elif method == 'map':
        return de_sitter_interval_map(N, D, KT2, fix_ends)
    else:
        assert False, 'Invalid method %s given to de_sitter_interval' % method

def de_sitter_interval_scatter(N, D, KT2, fix_ends=False):
    """ Scatter N points in a D dimensional interval in de Sitter spacetime

    This function uses the method described in Meyer1988 - a rejection method

    We scatter using a conformal factor of sigma=1+K/4(ds^2)

    """
    assert 0. < (KT2) < 4., 'KT^2 must be between 0 and 4 for this method'
    Z = np.empty(shape=(0, D))

    # rejection method
    # go in batches of size N
    while Z.shape[0] < N:
        R = minkowski_interval(N, D, fix_ends=fix_ends)
        R[:, 1:] -= 0.5 # fix back to 0 centre spatially
        M = (1. - (KT2 * 0.25))**(-D)  # maximum value
        m = np.random.rand(N) * M           # random assignments in that range
        S = (-1. * R[:,0]**2) + np.sum(R[:,1:]**2, axis=1) # proper time for each point
        sigma = (1. + (0.25 * KT2 * S))**(-D)
        Z = np.concatenate([Z, R[m < sigma]], axis=0)
    if fix_ends:
        Z[0,:] = 0.
        Z[1,:] = 0.
        Z[1,0] = 1.
    return Z[:N]

def de_sitter_interval_map(N, D, KT2, fix_ends=False):
    assert False, 'Not implemented yet'

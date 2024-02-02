"""
Models for causal sets
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

__all__ = ['CubeSpaceGraph', 'cube_space_interval']


class CubeSpaceGraph(AbstractDAG):
    def __init__(self, lp=1):
        super().__init__()
        self.lp = lp

    def _is_neighbor(self, x, y, weights):
        return (x < y).all(axis=1)

    def _calculate_weights(self, x, y):
        return dag.lp_distance(x, y, self.lp)
    

def cube_space_interval(N, D, fix_ends=True, sorted=True):
    """ Scatter N points in a D dimensional cube.
    Start point has coordinates: (0, 0, 0, ..., 0).
    End point has coordinates: (1, 1, 1, ..., 1).

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

    if fix_ends:
        a, b = np.zeros(D), np.ones(D)
        R[0], R[1] = a, b

    if sorted:
        R.sort(axis=0)

    return R

"""
Models for causal sets
"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)"])

import networkx as nx
import numpy as np


__all__ = ['cube_space_graph']

def cube_space_graph(N, D, p=1.0):
    """
    Create a cube space DAG


    Parameters
    ----------

    N - number of vertices in DAG
    D - dimension of box space
    p - probability with which allowed edges appear

    Notes
    -----

    In the cube space model, point a connects to point b iff a has smaller
    coordinates than b in every dimension.
    D=1 is a random DAG
    D=2 is equivalent to Minkowski space with D=2.
    """
    R = np.random.random((N, D))
    G = nx.DiGraph()
    edge_list = []
    for i in range(N):
        G.add_node(i, position=tuple(R[i]))
        for j in range(N):
            if (R[i] > R[j]).all():
                if p == 1. or p > np.random.random():
                    edge_list.append([j,i])
    G.add_edges_from(edge_list)
    return G

""" Embedding DAGs in Minkowski spacetime"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)"])

import networkx as nx
import numpy as np
import dagology as dag

__all__ = ['minkowski_embed']


def minkowski_embed(G, D, node_list=None):
    """ Embed a DAG in Minkowski spacetime

    We are using naive matrix methods currently but could upgrade in future


    """
    if not node_list:
        node_list = list(G.nodes())
    A = nx.adjacency_matrix(G, node_list).toarray()
    LP = dag.longest_path_matrix(A)
    ds2 = dag.naive_spacelike_matrix(LP)
    X = dag.mds(ds2, D, method='lorentzian')
    return X, node_list

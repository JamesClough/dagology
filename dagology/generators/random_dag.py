"""
Random DAG model, as in Karrer & Newman, 2009, Phys Rev E
"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)"])

import networkx as nx
import numpy as np
from random import randrange

import dagology as dag

__all__ = ['random_dag']

def random_dag(degree_sequence):
    """ Create a random DAG from a given degree sequence
    
    Parameters
    ----------
    
    degree_sequence - list of pairs of in, out degrees
    all edges go from earlier to later in this list
    
    Returns
    -------
    
    NetworkX DiGraph
    """
    G = nx.DiGraph()
    G.add_nodes_from(range(len(degree_sequence)))
    remaining_stubs = [] # list of forward pointing stubs
    for node, degrees in enumerate(degree_sequence):
        indegree, outdegree = degrees
        allowed_stubs = remaining_stubs[:]
        for x in range(indegree):
            if len(allowed_stubs) == 0:
                # raise networkx error
                assert False, 'Not a valid degree sequence'
            older_node = allowed_stubs.pop(randrange(len(allowed_stubs)))
            remaining_stubs.remove(older_node)
            # be careful about multiedges
            allowed_stubs = [x for x in allowed_stubs if x != older_node]
            G.add_edge(older_node, node)
        for x in range(outdegree):
            remaining_stubs.append(node)
    return G
    

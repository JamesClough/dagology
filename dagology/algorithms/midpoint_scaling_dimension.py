"""
Midpoint Scaling dimension estimator
"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)"])

import networkx as nx
import numpy as np
import dagology as dag

__all__ = ['mpsd']

def sub_interval_sizes(G, a, b, i):
    """
    Find the sizes of subintervals of an interval in a DAG
    
    Parameters
    ----------
    
    G : Networkx DiGraph
    a : Node in G
    b : Node in G, such that a < b (a precedes b in the graph)
    i : Node in G, such that a < i < b

    Returns
    -------
    
    2-tuple of the number of nodes in the [a,i] interval and the [i,b] interval
    """
    I_a = dag.interval(G, a, i)
    I_b = dag.interval(G, i, b)
    return (I_a.number_of_nodes(), I_b.number_of_nodes())


def mpsd(G):
    """
    Calculate the midpoint scaling dimension of a DAG

    Parameters
    ----------

    G : Networkx DiGraph
    """
    if G.number_of_edges() == 0:
        return 0.
    LP = nx.dag_longest_path(G)
    if len(LP) < 5:
        return 0.
    u, v = LP[0], LP[-1]

    I = dag.interval(G, u, v)

    # easy method - just check every item on the midpoint
    # hard method - start at the middle and check to see where value first drops
    # easy is implemented here for now
    max_N_min = 0
    max_intervals = None
    for i, w in enumerate(LP):
        intervals = sub_interval_sizes(I, u, v, w)
        N_min = min(intervals)
        if N_min > max_N_min:
            max_N_min = N_min
            max_intervals = intervals

    I_total = I.number_of_nodes()
    sub_I_total = sum(max_intervals) - 1.  # midpoint appears twice
    D = np.log2(I_total / sub_I_total)
    return (D + 1)

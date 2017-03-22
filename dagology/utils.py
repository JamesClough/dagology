""" Utility functions for DAG analysis"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)"])

import math
import networkx as nx
import numpy as np

__all__ = ['interval',
           'count_chains',
           'sphere_volume',
           'sphere_volume_analytic_cont']


def interval(G, a, b):
    """
    Calculate return the graph of the interval I[a,b] in G
    
    Parameters
    ----------
    
    G - NetworkX DiGraph (DAG)
    a - node in G
    b - node in G
    
    Returns
    -------
    
    I - NetworkX DiGraph - the interval [a,b]
    """
    a_dec = nx.descendants(G, a).union([a])
    b_anc = nx.ancestors(G, b).union([b])
    I_nodes = a_dec.intersection(b_anc)
    I = G.subgraph(I_nodes)
    return I

def count_chains(G, k):
    """
    Count the number of k-chains in G
    
    Parameters
    ----------
    
    G - NetworkX DiGraph (DAG)
    k - int - length of chains to count
    
    Returns
    -------
    
    C_k - number of chains of length k in G
    """
    sorted_nodes = nx.topological_sort(G)
    anc = {x:len(nx.descendants(G, x)) for x in sorted_nodes}
    for j in range(2, k):
        old_anc = anc
        anc = {x:int(np.sum([old_anc[y] for y in nx.descendants(G, x)])) for x in anc}        
    return np.sum(anc.values())

def sphere_volume(d, r=1.):
    """
    Volume bounded by d-Sphere
    d=0 is a line
    d=1 is a circle
    d=2 is a sphere etc.
    
    Parameters
    ----------
    
    d - int - dimension of sphere
    r - float - radius of sphere (default to unit sphere)"""
    assert d >= 0, "Sphere dimension must be positive"
    assert d == int(d), "Sphere dimension must be integer"
    assert r > 0, "Sphere radius must be positive"
    if d == 0:
        return 2. * r
    if d == 1:
        return (np.pi * r * r)
    else:
        return (sphere_volume(d-2, r) * r * r * 2 * np.pi / (d+1))

def sphere_volume_analytic_cont(d, r=1.):
    """
    Volume bounded by d-Sphere, for real d
    d=0 is a line
    d=1 is a circle
    d=2 is a sphere etc.
    
    Parameters
    ----------
    
    d - float - dimension of sphere
    r - float - radius of sphere (default to unit sphere)"""
    
    return  np.pi**((d+1.)/2.) * r**(d+1.) / (math.gamma((d+3.)/2.))
    
    
    
    


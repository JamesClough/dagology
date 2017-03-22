"""
Myrheim-Meyer Dimension estimator
"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)"])

import math       # For Gamma function
import networkx as nx
import numpy as np
import os
import dagology as dag

__all__ = ['mmd', 'mmd_formula', 'mmd_estimate']
D_MAX = 10

def mmd(G, k=2, already_tc=False):
    """
    Calculate the Myrheim-Meyer dimension of a DAG

    Parameters
    ----------

    G : Networkx DiGraph
    k : int
        Length of chains to count - default to 2
    """
    if G.number_of_edges() == 0:
        return 0

    if not already_tc:
        G = nx.transitive_closure(G)
    N = G.number_of_nodes()
    if k == 2:
        # this is a special case where we can use the inbuilt
        # number_of_edges function
        S = G.number_of_edges()
        f_D = float(S) / (N ** 2.0)
    else:
        S = dag.count_chains(G, k)
        f_D = float(S) / (N ** k)
        
    # lookup inverse of f_D(k) to find D estimate
    D = mmd_lookup(f_D, k)
    return D

def mmd_estimate(S, k, N):
    """ Estimate Myrheim-Meyer dimension from given number of k-chains
    
    Parameters
    ----------
    
    S : Number of k-chains
    k : Length of counted chains
    N : Number of vertices
    
    Returns
    -------
    D : Myrheim-Meyer dimension estimate
    """
    f_D = float(S) / (N ** k)
    D = mmd_lookup(f_D, k)
    return D

def mmd_lookup(f, k):
    D_range = np.arange(1, D_MAX, 0.01)
    for D in D_range:
        f_D = mmd_formula(D, k)
        if f_D < f:
            return D
    return D_MAX

def mmd_formula(D, k):
    """ Calculate Myrheim-Meyer f_D formula
    
    Parameters
    ----------
    
    D : Minkowski spacetime dimension (space + time)
    k : Length of chains being counted
    
    Returns
    -------
    
    S_k / N^k : Expected number of k-chains per N^k vertices in the interval
    For k=2 this is equal to the ordering fraction"""
    if D < 1:
        raise ValueError('D value %s less than 1' % D)
    if k==2:
        top = math.gamma(D + 1.) * math.gamma(D / 2.)
        bottom = 4. * math.gamma(D * 1.5)
    else:
        top = math.gamma(D / 2.) * math.gamma(D) * math.gamma(D+1.)**(k-1.)
        bottom = 2**(k-1) * k * math.gamma(0.5 * k * D) * math.gamma(0.5 * D * (k + 1))
        
    return top / bottom

def mmd_variance(D, k):
    """ Calculate variance of MMD estimate
    
    Currently only know this for k=2 but that's what we most commonly use"""
    assert k==2, 'Only known for k=2'
    gamma_arg = 1 + ((math.gamma(3.*D/2.) * 4)/(math.gamma(D/2.) * math.gamma(D+1)))
    var_c2 = 2. * mmd_formula(D, 3) * gamma_arg
    var_c2 += mmd_formula(D, 2)
    return var_c2


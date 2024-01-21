"""
De Sitter spacetime dimension and curvature estimator
"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)"])

__all__ = ['de_sitter_param_estimate']

import math
import networkx as nx
import numpy as np
import os
import dagology as dag
from scipy import optimize

def C_1(T, d, K, max_sum=10):
    """
    Calculate expected number of elements in unit density sprinkling into
    de Sitter space causal set of height T, dimension d+1, Gaussian curvature K

    Using formula from Meyer1988, 'Dimension of causal sets'
    """
    V = dag.sphere_volume_analytic_cont(d-1) / 2**d
    S = 0
    for i in range(max_sum):
        i_term = (K/4.)**i
        i_term *= T**((d+1) + 2*i)
        i_term *= G_1(d, i)
        S += i_term
    return S * V

def C_2(T, d, K, max_sum=10):
    """
    Calculate expected number of 2-chains in unit density sprinkling into
    de Sitter space causal set of height T, dimension d+1, Gaussian curvature K

    Using formula from Meyer1988, 'Dimension of causal sets'
    """
    V = (dag.sphere_volume_analytic_cont(d-1) / 2**d)**2
    S = 0
    for i_1 in range(max_sum):
        for i_2 in range(max_sum):
            i_term = (K / 4.)**(i_1 + i_2)
            i_term *= T**(2*(d + 1) + 2*(i_1 + i_2))
            i_term *= G_2(d, i_1, i_2)
            S += i_term
    return S * V

def C_3(T, d, K, max_sum=10):
    """
    Calculate expected number of 3-chains in unit density sprinkling into
    de Sitter space causal set of height T, dimension d+1, Gaussian curvature K

    Using formula from Meyer1988, 'Dimension of causal sets'
    """
    V = (dag.sphere_volume_analytic_cont(d-1) / 2**d)**3
    S = 0
    for i_1 in range(max_sum):
        for i_2 in range(max_sum):
            for i_3 in range(max_sum):
                i_term = (K / 4.)**(i_1 + i_2 + i_3)
                i_term *= T**(3*(d+1) + 2*(i_1 + i_2 + i_3))
                i_term *= G_3(d, i_1, i_2, i_3)
                S += i_term
    return S * V

def G_1(d, i):
    return 1./(2*i + d + 1)

def G_2(d, i_1, i_2):
    p = G_1(d, i_1) * 1./(2*(d + 1 + i_1 + i_2))
    p *= math.gamma(d + i_2 + 1)
    p /= math.gamma(i_2 + 1)
    p *= math.gamma(i_1 + i_2 + (d+3)/2.)
    p /= math.gamma(i_1 + i_2 + (d+1)*(3./2.))
    return p

def G_3(d, i_1, i_2, i_3):
    p = G_2(d, i_1, i_2) * 1./(3*(d+1) + 2*(i_1 + i_2 + i_3))
    p *= math.gamma(d + i_3 + 1)
    p /= math.gamma(i_3 + 1)
    p *= math.gamma(d + 2 + i_1 + i_2 + i_3)
    p /= math.gamma(2*d + 2 + i_1 + i_2 + i_3)
    return p

def equation_system(params, chains):
    T, d, K = params
    return (C_1(T, d, K) - chains[0],
            C_2(T, d, K) - chains[1],
            C_3(T, d, K) - chains[2])

def f(params, chains):
    return abs(sum(np.array(equation_system(params, chains))**2))

def find_T(K, N, d=1.):
    for T in np.arange(1., 100., 0.01):
        if C_1(T, d, K, max_sum=20) > N:
            break
    return T

def find_KT(KT2, N, d=1.):
    for T in np.arange(1., 100., 0.01):
        K = KT2 / (T*T)
        if dag.C_1(T, d, K, max_sum=20) > N:
            break
    return (K, T)

def de_sitter_param_estimate(chains, initial_guess=None, debug=False):
    """
    Estimate parameters of embedding de Sitter spacetime

    Parameters
    ----------

    chains - 3-tuple of numbers of k-chains for k=1,2,3
    initial_guess - guess of parameters (T, d, K) to seed optimiser
    debug - bool, print debugging statements

    Returns
    -------

    (T, d, K) - estimated embedding parameters

    """
    assert len(chains) == 3, 'Need to use N, C_2, C_3'

    if not initial_guess:
        # guess dimension assuming flat space as an initial guess
        D = dag.mmd_estimate(chains[1], 2, chains[0])
        d_guess = D - 1.
        if debug:
            print('d_guess = %s' % d_guess)
        initial_guess = [20., d_guess, 0.]

    min_opt = optimize.minimize(f, x0=initial_guess,
                                args=(chains,),
                                bounds=((0., 100.,), (0.01, 5.,), (-0.2, 0.2)))

    if debug:
        print(min_opt)
    T, d, K = min_opt['x']

    # check answer is reasonable
    if debug:
        diff_1 = (C_1(T, d, K) - chains[0]) / chains[0]
        diff_2 = (C_2(T, d, K) - chains[1]) / chains[1]
        diff_3 = (C_3(T, d, K) - chains[2]) / chains[2]
        if diff_1 > 0.1:
            print('C_1 out by %s' % diff_1)
        if diff_2 > 0.1:
            print('C_2 out by %s' % diff_2)
        if diff_3 > 0.1:
            print('C_3 out by %s' % diff_3)

    return T, d, K


if __name__ == "__main__":
    print(__doc__)

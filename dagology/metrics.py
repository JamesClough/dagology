""" Provide various metrics for calculating distances between two points
General form - metric(x, y) returns real number distance between points x and y

Note we can use these in
scipy.spatial.distance.pdist()

"""

#    Copyright (C) 2016 by
#    James Clough <james.clough91@gmail.com>
#    All rights reserved.
#    BSD license.

__author__ = "\n".join(["James Clough (james.clough91@gmail.com)", "Nicolas Kozak (nicolas.kozak5@gmail.com)"])

import numpy as np

def spherical(x_, y_):
    """Calculate distance on the surface of a d-sphere between points x and y

    We are using standard angular coordinates where the x[d-1] is in (0, 2*pi)
    and the other angles are in (0, pi)"""
    assert len(x_) == len(
        y_), 'ERROR - vectors in spherical metric have different lengths'
    if np.array_equal(x_, y_):
        return 0.
    if len(x_) == 1:
        # special case can be calculated more quickly
        return min(np.abs(x_[0] - y_[0]), (2. * np.pi - np.abs(x_[0] - y_[0])))
    x, y = angular_to_cartesian(x_), angular_to_cartesian(y_)
    cos_psi = np.dot(x, y)
    psi = np.arccos(cos_psi)
    return psi


def angular_to_cartesian(a):
    """Convert D angular spherical coordinates to D+1 cartesian - assume radius=1 """
    D = len(a)
    x = np.ones(D + 1)
    i = 0
    for i in range(D):
        x[i] *= np.cos(a[i])
        x[i + 1:] *= np.sin(a[i])
    return x


def cartesian_to_angular(a):
    """Convert d cartesian to d-1 angular spherical coordinates - assume radius=1 """
    D = len(a)
    x = np.zeros(D - 1)
    for i in range(D - 1):
        psi = np.arccos(a[i] / np.sqrt(np.sum((a[i:])**2)))
        x[i] = psi
    if a[-1] < 0:
        x[-1] = (2. * np.pi) - x[-1]
    return x


def hyperbolic(x, y, a=1.):
    """Calculate hyperbolic distances between coordinates in native representation """
    if np.array_equal(x, y):
        return 0.
    d_theta = spherical(x[1:], y[1:])
    cosh_ad = (np.cosh(a * x[0]) * np.cosh(a * y[0])) - \
        (np.sinh(a * x[0]) * np.sinh(a * y[0]) * np.cos(d_theta))
    d = np.arccosh(cosh_ad) / a
    return d * d


##########################################################################
# Lorentzian
##########################################################################

def minkowski(x, y, c=1.):
    """Calculate Minkowski separation between x and y using +--- convention
    c - speed of light - default to 1.
    y is allowed to be an array of size NxD
    The output is a 1D array of length N
    """
    x, y = np.asarray(x), np.asarray(y)
    if y.ndim == 1:
        y = y.reshape(1, -1)
    dx = y[:, 1:] - x[1:]
    return (c*(y[:, 0] - x[0]))**2 - (dx**2).sum(axis=1)


def minkowski_periodic(x, y, period, c=1.):
    """Calculate Minkowski separation between x and y using -++...+ convention
       Periodic boundary conditions in spatial coordinates are given in L
       If len(L) < D-1 then assume no boundary on other spatial dimensions"""
    D = len(x)
    assert len(y) == D, 'ERROR - vectors in minkowski have different lengths'
    if np.array_equal(x, y):
        return 0.
        
    while len(period) < (D - 1):
        period.append(None)
    dt = x[0] - y[0]
    dt2 = dt * dt
    ds2 = -1 * dt2
    for d in range(1, D):
        period_d = period[d - 1]  # the first dimension in the coords is time so exclude it
        if period_d:
            dx2 = min((x[d] - y[d])**2, (x[d] - y[d] + period_d)
                      ** 2, (x[d] - y[d] - period_d)**2)
        else:
            dx2 = (x[d] - y[d])**2
        ds2 += dx2
    return ds2


def de_sitter(x, y):
    """ Calculate de Sitter separation between x and y in conformal coordinates"""
    assert len(x) == len(
        y), 'ERROR - vectors in de Sitter metric have different lengths'
    if np.array_equal(x, y):
        return 0.
    dt = x[0] - y[0]
    dt2 = dt * dt
    dx = spherical(x[1:], y[1:])
    dx2 = dx * dx
    return (dx2 - dt2)


##########################################################################
# Lp norm
##########################################################################

def lp_distance(x, y, p):
    x, y = np.asarray(x), np.asarray(y)
    if y.ndim == 1:
        y = y.reshape(1, -1)
    return (abs(y - x)**(p)).sum(axis=1)
    
if __name__ == "__main__":
    print(__doc__)

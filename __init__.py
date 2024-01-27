""" Dagology

Dagology is a package of algorithms for analysing Directed Acyclic Graphs,
or DAGs by considering their causal structure, and the correspondence
with Lorentzian Geometry.

Find theoretical details in these papers which use these techniques:

What is the dimension of citation space?
http://www.sciencedirect.com/science/article/pii/S037843711501081X

Transitive reduction of citation networks
https://comnet.oxfordjournals.org/content/3/2/189.full

Lorentzian Multidimensional Scaling
http://arxiv.org/pdf/1602.03103.pdf
"""

import sys
if sys.version_info[:2] < (3, 6):
    m = "Python version 3.6 or later is required (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from dagology import *

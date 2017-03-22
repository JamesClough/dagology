from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true
from nose.tools import assert_almost_equal

import networkx as nx
import numpy as np
import dagology as dag
from scipy.misc import comb

class TestCountChains(object):
    """ Unit tests for chain counting function"""
    
    def test_unconnected(self):
        G = nx.DiGraph()
        G.add_nodes_from([1, 2, 3])
        S = dag.count_chains(G, 2)
        assert_equal(S, 0)
        S = dag.count_chains(G, 3)
        assert_equal(S, 0)

    def test_line(self):
        G = nx.path_graph(10, create_using=nx.DiGraph())
        for k in range(2, 11):
            S = dag.count_chains(G, k)
            assert_equal(S, comb(10, k))
            
    def test_complete(self):
        G = nx.path_graph(10, create_using=nx.DiGraph())
        G = nx.transitive_closure(G)
        for k in range(2, 11):
            S = dag.count_chains(G, k)
            assert_equal(S, comb(10, k))
            
    def test_example(self):
        G = nx.DiGraph()
        G.add_edges_from([[1,2], [1,3], [2,4], [3,5], [4,6], [5,6]])
        assert_equal(dag.count_chains(G, 2), 11)
        assert_equal(dag.count_chains(G, 3), 8)
        assert_equal(dag.count_chains(G, 4), 2)
        assert_equal(dag.count_chains(G, 5), 0)

class TestInterval(object):
    """ Unit tests for interval function"""

    def test_link(self):
        G = nx.DiGraph()
        G.add_edge(1, 2)
        I = dag.interval(G, 1, 2)
        assert_equal(list(G.nodes()), list(I.nodes()))
        assert_equal(list(G.edges()), list(I.edges()))

    def test_unconnected(self):
        G = nx.DiGraph()
        G.add_nodes_from([1, 2, 3])
        I = dag.interval(G, 1, 2)
        assert_equal(I.number_of_nodes(), 0)

    def test_example(self):
        G = nx.DiGraph()
        G.add_edges_from([[1, 2], [1, 3], [1, 4], [2, 5], [3, 5], [5, 6]])
        I = dag.interval(G, 1, 5)
        assert_equal(sorted(list(I.nodes())), [1, 2, 3, 5])

class TestSphereVolume(object):
    """ Unit tests for sphere voume function"""
    
    def test_d0(self):
        assert_equal(dag.sphere_volume(0), 2.)
        assert_equal(dag.sphere_volume(0, 3), 6.)
        
    def test_d1(self):
        assert_almost_equal(dag.sphere_volume(1), np.pi)
        assert_almost_equal(dag.sphere_volume(1, 10), 100 * np.pi)

    def test_d2(self):
        assert_almost_equal(dag.sphere_volume(2), (4./3.) * np.pi)
        assert_almost_equal(dag.sphere_volume(2, 10), (4./3.) * 1000 * np.pi)
        
class TestSphereVolumeAnalyticCont(object):
    """ Unit tests for sphere voume function when analytically continued"""
    
    def test_equality_real_d(self):
        for d in [1, 2, 3, 4, 5]:
            for r in [0.1, 1., 2., 10.,]:
                assert_almost_equal(dag.sphere_volume(d, r), 
                                     dag.sphere_volume_analytic_cont(d, r))

               
        

    

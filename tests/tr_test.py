""" Test transitive reduction"""
import dagology as dag
import unittest
import numpy as np

class TestTransitiveReduction(unittest.TestCase):
    def test_chain(self):
        """ Create fully connected N node graph and check transitive reduction is a chain"""
        N = 10
        A = np.zeros((N,N))
        for i in range(N):
            for j in range(i):
                A[i,j] = 1
        A_tr = dag.transitive_reduction(A)
        for i in range(N):
            for j in range(N):
                if i==j+1:
                    self.assertEqual(A_tr[i,j], 1)
                else:
                    self.assertNotEqual(A_tr[i,j], 1)

if __name__ == "__main__":
    unittest.main()
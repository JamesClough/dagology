""" Test transitive completion"""
import numpy as np
import dagology as dag
import unittest

class TestTransitiveCompletion(unittest.TestCase):
    def test_chain(self):
        """ Create chain N node graph and check transitive completion"""
        N = 10
        A = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                if i == j+1:
                    A[i,j] = 1
        A_tc = dag.transitive_completion(A)
        for i in range(N):
            for j in range(N):
                if i>j:
                    self.assertEqual(A_tc[i,j], 1)
                else:
                    self.assertNotEqual(A_tc[i,j], 1)

if __name__ == "__main__":
    unittest.main()

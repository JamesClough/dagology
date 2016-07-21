import unittest
import dagology as dag
import numpy as np

class TestLongestPath(unittest.TestCase):
    """ Test Longest path matrix"""
    def test_longest_path_chain(self):
        N = 10
        A = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                if i == j+1:
                    A[i,j] = 1
        self.assertEqual(np.max(dag.longest_path_matrix(A)), N-1)

	def test_longest_path_chain(self):
		N = 10
        A = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                if i == j+1:
                    A[i,j] = 1
        self.assertEqual(np.max(dag.longest_path_matrix(A)), N-1)

class TestNaiveSpacelikeSeparation(unittest.TestCase):
    """ Test Naive spacelike distance matrix
    TODO: WRITE THESE TESTS"""
    pass

if __name__ == "__main__":
    unittest.main()

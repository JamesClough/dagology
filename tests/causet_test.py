import unittest
import dagology as dag
import numpy as np

class TestScatter(unittest.TestCase):
    def test_shape(self):
        N = 100
        for D in range(1,10):
            for opt in [True, False]:
                A, R = dag.minkowski_causet(N, D, fix_ends=opt)
                # check coordinates
                self.assertEqual(R.shape, (N,D))
                self.assertTrue(np.max(R) <= 1.)
                self.assertTrue(np.min(R) >= -1.)
                # check adjacency matrix
                self.assertEqual(A.shape, (N,N))
                self.assertEqual(np.max(A), 1)
                self.assertEqual(np.min(A), 0)


if __name__ == "__main__":
    unittest.main()
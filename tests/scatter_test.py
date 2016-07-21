import dagology as dag
import numpy as np
import unittest

class TestScatter(unittest.TestCase):
    def test_spherical(self):
        N = 100
        R = dag.sphere_surface_cartesian(N, 5)
        X = dag.cartesian_to_angular(R[0])
        C = dag.angular_to_cartesian(X)
    	self.assertTrue(np.max(np.abs(R[0] - C)) < 0.0001)

if __name__ == "__main__":
    unittest.main()

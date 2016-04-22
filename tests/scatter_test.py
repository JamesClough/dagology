import dagology as dag
import numpy as np
import matplotlib.pyplot as plt

def test_spherical_scatter():
    R = dag.sphere_surface_cartesian(1, 5)
    X = dag.cartesian_to_spherical(R[0])
    C = dag.spherical_to_cartesian(X)
    assert np.max(np.abs(R[0] - C)) < 0.0001
    
def de_sitter_scatter():
    R = dag.de_sitter_interval_scatter(10000, 2, 0.1, 1.3)
    plt.scatter(np.sin(R[:2, 1]), R[:2, 0], color='r')
    plt.scatter(np.sin(R[2:, 1]), R[2:, 0], color='b')
    plt.show()
    
if __name__ == "__main__":
    de_sitter_scatter()

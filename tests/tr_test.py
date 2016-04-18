""" Test transitive reduction"""

import numpy as np
import dagology as dag

def test_tr():
    A = np.zeros((10, 10))
    N, _ = A.shape
    for i in range(N):
        for j in range(N):
            if i > j:
                A[i,j] = 1
    A_tr = dag.transitive_reduction(A)
    print A
    print A_tr
    
def main():
    test_tr()
    
if __name__ == "__main__":
    main()

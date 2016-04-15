""" Test transitive completion"""
import numpy as np
import dagology as dag

def test_tc():
    A = np.zeros((10, 10))
    N, _ = A.shape
    for i in range(N):
        for j in range(N):
            if i == j+1:
                A[i,j] = 1
    A_tc = dag.transitive_completion(A)
    print A
    print A_tc
    
def main():
    test_tc()
    
if __name__ == "__main__":
    main()

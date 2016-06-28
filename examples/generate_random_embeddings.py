import dagology as dag
import example_utils as util
import pickle 
import networkx as nx

def random_dag_data(N, k, D_hat, r, base_filepath):
    for i in range(r):
        print 'Creating random dag %s of %s' % (i+1,r)
        A = dag.random_dag(N, k)
        LP = dag.longest_path_matrix(A)
        X = dag.naive_spacelike_matrix(LP)
        R = dag.MDS(X, D_hat, method='lorentzian')
        R = util.flip_time_axis(R, A)
        fp = base_filepath + ('_%s_R.pkl' % i)
        with open(fp, 'wb') as f:
            pickle.dump(R, f)
        fp = base_filepath + ('_%s_A.pkl' % i)
        with open(fp, 'wb') as f:
            pickle.dump(A, f)

def causal_set_data(N, D, D_embed, r, base_filepath):
    for i in range(r):
        print 'Creating causal set %s of %s' % (i+1,r)
        A, R = dag.minkowski_causet(N, D)
        LP = dag.longest_path_matrix(A)
        SP = dag.naive_spacelike_matrix(LP)
        X = dag.MDS(SP, D, method='lorentzian')
        X = util.flip_time_axis(X, A)
        if D==2:
            X = util.flip_space_axis(X, R, 1)
        R -= 0.5
        X = util.normalise_coords(X)
        fp = base_filepath + ('_%s_R.pkl' % i)
        with open(fp, 'wb') as f:
            pickle.dump(R, f)
        fp = base_filepath + ('_%s_A.pkl' % i)
        with open(fp, 'wb') as f:
            pickle.dump(A, f)
        fp = base_filepath + ('_%s_X.pkl' % i)
        with open(fp, 'wb') as f:
            pickle.dump(X, f)
    
if __name__ == "__main__":
    'https://docs.python.org/2.7/library/argparse.html' # for later
    import sys
    try:
        dag_type = sys.argv[1]
    except:
        dag_type = 'r'
        
    if dag_type == 'r':
        N = int(sys.argv[2])
        k = int(sys.argv[3])
        r = int(sys.argv[4])
        base_filepath = '../data/random_dag/N%s_k%s' % (N, k)
        random_dag_data(N, k, 10, r, base_filepath)
        
    if dag_type == 'c':
         N = int(sys.argv[2])
         D = int(sys.argv[3])
         r = int(sys.argv[4])
         base_filepath = '../data/causal_sets/N%s_D%s' % (N, D)
         causal_set_data(N, D, 10, r, base_filepath)
         

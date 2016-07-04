import dagology as dag
import example_utils as util
import pickle 
import numpy as np

def embed_citnet(name, N, D):
    filepath = './datasets/%s.txt' % name
    edges, num_vertices = util.load_txt_dataset(filepath) 
    top_N_edges, top_N_nodes = util.restrict_top_N(edges, N)
    A = np.zeros((N, N))
    K = np.zeros(N)
    for a, b in top_N_edges:
        A[a,b] = 1 
        K[b] += 1
    #A_tr = dag.transitive_reduction(A)
    LP = dag.longest_path_matrix(A)
    max_sep = int(np.max(LP)/2)
    SP = dag.naive_spacelike_matrix(LP, max_sep)
    R = dag.MDS(SP, D, method='lorentzian')
    R = util.flip_time_axis(R, A)
    A_filepath = './results/%s_%s_A.pkl' % (name, N)
    X_filepath = './results/%s_%s_X.pkl' % (name, N)
    K_filepath = './results/%s_%s_K.pkl' % (name, N)
    nodes_filepath = './results/%s_%s_nodes.pkl' % (name, N)
    with open(A_filepath, 'wb') as f:
        pickle.dump(A, f)
    with open(X_filepath, 'wb') as f:
        pickle.dump(R, f)
    with open(K_filepath, 'wb') as f:
        pickle.dump(K, f)
    with open(nodes_filepath, 'wb') as f:
        pickle.dump(top_N_nodes, f)

def load_citnet(name, N):
    A_filepath = './results/%s_%s_A.pkl' % (name, N)
    X_filepath = './results/%s_%s_X.pkl' % (name, N)
    K_filepath = './results/%s_%s_K.pkl' % (name, N)
    nodes_filepath = './results/%s_%s_nodes.pkl' % (name, N)
    with open(A_filepath, 'rb') as f:
        A = pickle.load(f)
    with open(X_filepath, 'rb') as f:
        X = pickle.load(f)
    with open(K_filepath, 'rb') as f:
        K = pickle.load(f)
    with open(nodes_filepath, 'rb') as f:
        top_N_nodes = pickle.load(f)
    return {'A':A,
            'X':X,
            'K':K,
            'nodes':top_N_nodes}
    
        
if __name__ == "__main__":
    import sys
    name, N = sys.argv[1:3]
    N = int(N)
    try:
        D = sys.argv[3]
        D = int(D)
    except:
        D = 10
    embed_citnet(name, N, D)


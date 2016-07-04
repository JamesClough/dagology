# useful definitions and functions for dagology examples
import numpy as np

def load_txt_dataset(filepath):
    edges = []
    from sets import Set
    vertices = Set([])
    with open(filepath, 'r') as f:
        for line in f:
            a, b = line.split()
            edges.append((int(a),int(b)))
            vertices.add(int(a))
            vertices.add(int(b))
    vertices = list(vertices)
    vertices.sort()
    vertex_dict = {v:i for i, v in enumerate(vertices)}
    index_edges = [(vertex_dict[a], vertex_dict[b]) for a, b in edges]
    return index_edges, len(vertices)

# restrict to N nodes with largest degree
def restrict_top_N(edges, N):
    from collections import defaultdict
    degree = defaultdict(int)
    for a,b in edges:
            #degree[a] += 1
            degree[b] += 1
    degree_count = [(degree[x], x) for x in degree.keys()]
    degree_count.sort()
    degree_count.reverse()
    top_N = [x[1] for x in degree_count[:N]]
    top_N.sort()
    top_N_edges = []
    for a,b in edges:
        if a in top_N:
            if b in top_N:
                top_N_edges.append((top_N.index(a),top_N.index(b)))
    return top_N_edges, top_N

def flip_time_axis(R, A):
    N, _ = A.shape
    right_edges = 0
    wrong_edges = 0
    for i in range(N):
        for j in range(N):
            if A[i,j]:
                if R[i,0]>R[j,0]:
                    right_edges += 1
                else:
                    wrong_edges += 1
    if wrong_edges > right_edges:
        R[:,0] *= -1.
    return R

def flip_space_axis(X, R, d=1):
    """ Flip d-th spatial axis of X to make it more similar to R"""
    same = 0
    N = R.shape[0]
    for i in range(N):
        same += X[i,d] * R[i,d]
    if same < 0:
        X[:,d] *= -1.
    return X

def normalise_coords(X):
    """ Try to fit coords into a box around unit size"""
    height = np.max(X[:,0]) - np.min(X[:,0])
    X /= height
    return X
    
def average_deviation(R, X):
    """ Calculate average deviation of embedded point from where it should be"""
    N = R.shape[0]
    sq_dev = []
    for i in range(N):
        dev = R[i,:] - X[i, :]
        sq_dev.append(np.sum(dev**2))
    return np.sqrt(np.mean(sq_dev))

def measure_embedding(R, A, c=1.):
    # A needs to be transitively complete here
    true_pos = 0.
    true_neg = 0.
    false_pos = 0.
    false_neg = 0.
    N, _ = A.shape
    for i in range(N):
        for j in range(i):
            dR2 = (R[i,:] - R[j,:])**2
            ds2 = np.sum(dR2[1:])  - (c * dR2[0])
            if (ds2<0) and (R[i,0]>R[j,0]):
                lightcone = True
            else:
                lightcone = False
            if A[i,j] == 1:
                if lightcone:
                    true_pos += 1
                else:
                    false_neg += 1
            else:
                if lightcone:
                    false_pos += 1
                else:
                    true_neg += 1
    return true_pos, true_neg, false_pos, false_neg

def receiver_operator_curve(R, A, c=1.):
    true_pos, true_neg, false_pos, false_neg = measure_embedding(R, A, c)
    if true_pos == 0:
        sensitivity = 0
    else:
        sensitivity = true_pos / (true_pos + false_neg)
    if true_neg == 0:
        specificity = 0
    else:
        specificity = true_neg / (true_neg + false_pos)
    return sensitivity, specificity

def create_ROC_data(X, A):
    c_range = np.logspace(-2., 2., 40)
    x, y = [0.], [0.]
    for c in c_range:
        sens, spec = receiver_operator_curve(X, A, c)
        x.append(1.-spec)
        y.append(sens)
    x.append(1.)
    y.append(1.)
    return x, y
    
    
# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)    

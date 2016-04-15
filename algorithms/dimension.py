""" Dimension measures for DAGs

myrheim_meyer(A, k) - Myrheim-Meyer dimension for k chains

midpoint_scaling(A, i, j) - Midpoint-scaling dimension in I[i,j] interval
"""

def get_mmd_data(k, filepath=None):
    if filepath == None:
        filepath = '../cache/MMD_data/f_of_d_%s.txt' % k
    try:
        f_dict = {}
        with open(filepath) as f:
            for line in f.readlines():
                d, f_d = line.split()
                f_dict[float(f_d)] = float(d)
    except:
        print """myrheim_meyer Couldn't find MMD data for k=%s, in file %s so 
               creating new data now"""
        f_dict = create_mmd_data(k, filepath)
        
def create_mmd_data(k, filepath):
    pass

def ordering_fraction(A, k):
    """ Calculate the k-path ordering fraction"""
    A_k = A.copy(deep=True)
    for i in range(k-2):
        A_k = np.dot(A_k, A)
    N, _ = A.shape
    of = np.sum(A_k) / N**k
    return of

def myrheim_meyer(A, k=2, cached_f_dict=None, filepath=None):
    """ Calculate Myrheim-Meyer dimension counting k-chains
    
    IMPORTANT - this function assumes the adjacency matrix is transitively
                complete!"""
    
    # try to load MMD data - JC: will want to make this find the file later
    if cached_f_dict == None:
        f_dict = get_mmd_data(k, filepath)
    else:
        f_dict = cached_f_dict

    of = ordering_fraction(A, k)
    
    best_d = None
    min_diff = np.inf
    for f_d in f_dict:
        diff = abs(f_d - of)
        if diff < min_diff:
            min_diff = diff
            best_d = f_dict[f_d]
    return best_d
        
def lp(A, i=1, j=0):
    """ Calculate length of longest path from vertex i to j
    
    A -- adjacency matrix
    i -- index of starting vertex
    j -- index of final vertex
    Assumes we only have the [i,j] interval"""
    n = 1
    P = A.copy()            # this will be matrix of longest path lengths from i
    A_n = A.copy()
    while np.sum(A_n) > 0:
        C  = A_n.copy()     # connection matrix
        C[C > 1] = 1
        C *= n
        P = np.maximum(P, C)
        A_n = np.dot(A_n, A)
        n += 1
    return (n-1)

def interval(A, i=1, j=0):
    """ Return adjacency matrix of interval between vertex i and j"""
    i_set = Set([i])
    j_set = Set([j])
    A_n = A.copy()
    N, _ = A.shape
    while np.sum(A_n) > 0:
        for k in range(N):
            if A[i,k] > 0:
                i_set.add(k)
            if A[k,j] > 0:
                j_set.add(k)
        A_n = np.dot(A_n, A)
    i_set.intersection_update(j_set)
    A_int = A.copy()
    for k in range(N):
        if k not in i_set:
            A_int[k,:] = 0
            A_int[:,k] = 0
    return A_int
  
def get_midpoint(A, i, j):
    """ Find midpoint vertex
    
    Find midpoint vertex defined as the vertex k where the minimum of 
    the size of intervals I[i,k] and I[j,k] is maximised
    """
    N, _ = A.shape
    max_min_V = 0
    k_max = None
    for k in range(N):
        # select midpoint
        N1, N2 = interval(A, i, k), interval(A, k, j)
        V1, V2 = np.sum(N1[i, :]), np.sum(N2[:, j])
        min_V = min(V1, V2)
        if min_V > max_min_V:
            max_min_V = min_V
            k_max = k
    return k_max
    
def midpoint_scaling(A, i=1, j=0):
    """ Calculate Midpoint-scaling dimension"""
    k = get_midpoint(A)
    N, _ = A.shape
    A1, A2 = interval(A, i, k), interval(A, k, j)
    N1, N2 = np.sum(N1[i, :]), np.sum(N2[:, j])
    N_mean = np.mean(N1, N2)
    D = (np.log(N_mean) - np.log(N)) / np.log(2)
    return D
    

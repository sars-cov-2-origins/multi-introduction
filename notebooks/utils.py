import pandas as pd
import networkx as nx
import numpy as np
from networkx.algorithms import components

# issues with some sites
to_remove = [x for x in range(29700,29720)]

# read FASTA stream and return (ID,seq) dictionary
def readFASTA(path):
    stream = open(path)
    seqs = {}
    name = None
    seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq.upper()
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq.upper()
    return seqs


# compare 2 sequences and return the mutations (e.g., C8782T)
def compare_seqs(ref, seq):
    mut_indices = []
    for index, nt in enumerate(seq):
        if (ref[index] != seq[index]) and ref[index] in 'ACGT' and seq[index] in 'ACGT':
            mut_indices.append(ref[index] + str(index + 1) + seq[index])
    return mut_indices


# get lineage of Hu-1 aligned sequence
def get_lineage(seq):
    if seq[8782-1] == 'C' and seq[28144-1] == 'T':
        lineage = 'B'
    elif seq[8782-1] == 'T' and seq[28144-1] == 'C':
        lineage = 'A'
    elif seq[8782-1] == 'C' and seq[28144-1] == 'C':
        lineage = 'CC'
    elif seq[8782-1] == 'T' and seq[28144-1] == 'T':
        lineage = 'TT'
    else:
        lineage = 'unknown'
    return lineage


# get components of an alignment based on lineages and a reference sequence
def get_components(lin1, lin2, fasta, ref_seq, mask_pos = None, max_muts = None):
    # lin1 and lin2 must be chosen from {'linA', 'linB', 'CC', 'TT', 'unknown'}
    
    # create list of mutations of interest
    intermediate_mut_indices = []
    for key in fasta:
        seq = fasta[key]
        
        if lin1 == 'CC' or lin2 == 'CC':
            if seq[28144-1] == 'C' and seq[8782-1] == 'C':
                intermediate_mut_indices += compare_seqs(ref=ref_seq, seq=seq)
        elif lin1 == 'TT' or lin2 == 'TT':
            if seq[28144-1] == 'T' and seq[8782-1] == 'T':
                intermediate_mut_indices += compare_seqs(ref=ref_seq, seq=seq)
        else:
            if (seq[28144-1] == 'C' and seq[8782-1] == 'T') or (seq[28144-1] == 'T' and seq[8782-1] == 'C'):
                intermediate_mut_indices += compare_seqs(ref=ref_seq, seq=seq)

    intermediate_mut_indices = list(set(intermediate_mut_indices))
    # don't want these 2 mutations; can just call them lineage A/B/CC/TT as needed
    if 'T28144C' in intermediate_mut_indices:
        intermediate_mut_indices.remove('T28144C')
    if 'C8782T' in intermediate_mut_indices:
        intermediate_mut_indices.remove('C8782T')


    intermediate_homoplasy_d = {mut:{'linA':[], 'linB':[], 'CC':[], 'TT':[], 'unknown':[]} for mut in intermediate_mut_indices}
    seq_homoplasy_d = {}

    for key in fasta:
        seq = fasta[key]

        if (seq[8782-1] == 'C' and seq[28144-1] == 'T'):
            lineage = 'linB'
        elif (seq[8782-1] == 'T' and seq[28144-1] == 'C'):
            lineage = 'linA'
        elif (seq[8782-1] == 'C' and seq[28144-1] == 'C'):
            lineage = 'CC'
        elif (seq[8782-1] == 'T' and seq[28144-1] == 'T'):
            lineage = 'TT'
        else:
            lineage = 'unknown'

        for mut in intermediate_mut_indices:
            mut_index = int(mut[1:-1])
            sub = mut[-1]

            if mask_pos:
                if mut_index in mask_pos:
                    continue

            if seq[mut_index-1] == sub:
                intermediate_homoplasy_d[mut][lineage].append(key)

                if lineage == lin1 or lineage == lin2:
                    if key not in seq_homoplasy_d:
                        seq_homoplasy_d[key] = {'lineage':lineage, 'muts':[mut]}
                    else:
                        seq_homoplasy_d[key]['muts'].append(mut)

    
    # get list of mutations that are not in both lineages to avoid them later
    muts_to_remove = []
    for mut in intermediate_homoplasy_d:
        if (len(intermediate_homoplasy_d[mut][lin1]) > 0) and (len(intermediate_homoplasy_d[mut][lin2]) > 0):
            continue
        else:
            muts_to_remove.append(mut)

        if int(mut[1:-1]) in to_remove:
            muts_to_remove.append(mut)
    
    # make graph of sequences and mutations
    g = nx.Graph()
    for key in seq_homoplasy_d:
        if max_muts:
            if len(seq_homoplasy_d[key]['muts']) > max_muts:
                continue
        for mut in seq_homoplasy_d[key]['muts']:
            if mut in muts_to_remove:
                continue
            g.add_edge(key, mut)

    # make dictionary of components
    component = 1
    component_dict = {}
    for n in components.connected_components(g):
        skip = False
        comp_seqs = []
        comp_muts = []
        for _ in n:
            if 'EPI' in _:
                comp_seqs.append(_)
            else:
                comp_muts.append(_)
        if len(comp_seqs) == 1:
            continue
        component_dict[component] = {'seqs':comp_seqs, 'muts':comp_muts}
        component += 1   
        
        
    # pool sites together
    comp_sites = {}
    comp_df = pd.DataFrame()
    comp_data = []
    for comp in component_dict:
        comp_sites[comp] = []
        temp_comp_data = []
        for key in component_dict[comp]['seqs']:
            comp_sites[comp] += [x for x in seq_homoplasy_d[key]['muts'] if x not in muts_to_remove]
            temp_comp_data = temp_comp_data + [[key.split('|')[1], seq_homoplasy_d[key]['lineage']] + [x for x in seq_homoplasy_d[key]['muts'] if x not in muts_to_remove]]
            temp_comp_df = pd.DataFrame(temp_comp_data)
            temp_comp_df = temp_comp_df.sort_values(by=temp_comp_df.shape[1]-1).sort_values(by=1)
        comp_sites[comp] = set(comp_sites[comp])
    return component_dict, comp_sites, seq_homoplasy_d, muts_to_remove, intermediate_mut_indices


def random_seq(k, gtr_probs):
    from random import uniform
    '''
    This function generates a random sequence of length k using GTR stationary probabilities
    :param k: The length of the output sequence
    :param gtr_probs: The GTR stationary probabilities as a list [prob_A, prob_C, prob_G, prob_T]
    :return: A random string of length k generated using gtr_probs
    
    prob_A, prob_C, prob_G, prob_T = gtr_probs # You can use these if it's more convenient
    '''
    # TODO Your code here
    # set up variables/parameters
    letters = 'ACGT'
    gtr_scale = [sum(gtr_probs[:i]) for i in range(1, len(gtr_probs)+1)]
    seq = ''
    
    # loop through length of k and get a new probability each time
    for _ in range(k):
        prob = uniform(0, 1)
        index = 0
        while gtr_scale[index] < prob:
            index += 1
        seq += letters[index]
    return seq
    

def evolve(tree, root_seq, gtr_probs, gtr_rates):
    import treeswift
    import numpy as np
    from scipy import linalg
    import scipy as scp
    from random import uniform
    '''
    This function simulates the evolution of a root sequence down a given tree under the GTR model
    :param tree: A TreeSwift Tree object representing the phylogenetic tree
    :param root_seq: The root sequence to evolve down the tree
    :param gtr_probs: The GTR stationary probabilities as a list [prob_A, prob_C, prob_G, prob_T]
    :param gtr_rates: The GTR transition rates as a list [rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG]
    :return: A dictionary where keys are the labels of the given tree and values are evolved sequences
    '''
    letters = "ACGT"
    probs_matrices = []
    # set probs and rates
    rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG = gtr_rates
    prob_A, prob_C, prob_G, prob_T = gtr_probs
    seqs = dict() # This will be your output (keys = leaf labels (str) and values = sequences (str))
    
    # matrices from probs and rates
    pi_mat = np.identity(4) * [prob_A, prob_C, prob_G, prob_T]
    trans_mat = [[0] + [rate_AC, rate_AG, rate_AT], 
                 [rate_AC] + [0] + [rate_CG, rate_CT], 
                 [rate_AG, rate_CG] + [0] + [rate_GT],
                 [rate_AT, rate_CT, rate_GT] + [0]]
    
    # make the proper R matrix
    R = np.matmul(trans_mat, pi_mat) 
    R = [-sum(R[i1]) if i1 == i2 else R[i1][i2] for i1 in range(4) for i2 in range(4)] # creating the v values, except into a 1d array
    R = np.reshape(R, (4,4)) # reshaping the array
    
    # get the relevant D matrices and construct X
    D_pos12 = scp.linalg.fractional_matrix_power(pi_mat, 1/2)
    D_neg12 = scp.linalg.fractional_matrix_power(pi_mat, -1/2)
    X = np.linalg.multi_dot([D_pos12, R, D_neg12])
    # sanity check
    R1 = np.linalg.multi_dot([D_neg12, X, D_pos12])

    # eigen decomp of X
    u, s, vh = np.linalg.svd(X)
    smat = np.diag(s)
    ut = -u.T

    # getting D
    d = 0
    for i in range(4):
        d -= pi_mat[i][i]*R[i][i]
        
    # scale R
    R_scaled = R / d

    # need e^t*R_scaled[i,j] for each component now = matrix exponential
    # R_t = scp.linalg.expm(t*R_scaled)
    
    allSeqs = dict()
    # traverse the leaves
    for n in tree.traverse_preorder():
        # if root, assign seq
        if n.is_root():
            allSeqs[n] = root_seq
            continue
        
        # set edge length as distance and use that to make the prob matrix
        t = n.edge_length
        probs = scp.linalg.expm(t*R_scaled)
        # make a scale for each letter probability so i don't have to recompute this
        probs_scales = [[sum(letterProbs[:i]) for i in range(1, len(letterProbs)+1)] for letterProbs in probs]
        probs_matrices.append((n.label, n.edge_length, probs_scales))
        parentSeq = allSeqs[n.parent]
        childSeq = ""
        for nt in parentSeq:
            index1 = letters.index(nt)
            index2 = 0
            prob = uniform(0, 1)
            while probs_scales[index1][index2] < prob:
                index2 += 1
            childSeq += letters[index2]
        allSeqs[n] = childSeq
        
    for seq in allSeqs:
        if seq.is_leaf():
            seqs[seq.label] = allSeqs[seq]
        
    return seqs, probs_matrices


def calc_min_interval(x, alpha):
    """Internal method to determine the minimum interval of a given width
    Assumes that x is sorted numpy array.
    """

    n = len(x)
    cred_mass = 1.0-alpha

    interval_idx_inc = int(np.floor(cred_mass*n))
    n_intervals = n - interval_idx_inc
    interval_width = x[interval_idx_inc:] - x[:n_intervals]

    if len(interval_width) == 0:
        raise ValueError('Too few elements for interval calculation')

    min_idx = np.argmin(interval_width)
    hdi_min = x[min_idx]
    hdi_max = x[min_idx+interval_idx_inc]
    return hdi_min, hdi_max


def hpd_single(x, alpha):
    """Calculate highest posterior density (HPD) of array for given alpha. 
    The HPD is the minimum width Bayesian credible interval (BCI).
    :Arguments:
        x : Numpy array
        An array containing MCMC samples
        alpha : float
        Desired probability of type I error (defaults to 0.05)
    """

    # Make a copy of trace
    x = x.copy()
    # For multivariate node
    if x.ndim > 1:
        # Transpose first, then sort
        tx = np.transpose(x, list(range(x.ndim))[1:]+[0])
        dims = np.shape(tx)
        # Container list for intervals
        intervals = np.resize(0.0, dims[:-1]+(2,))

        for index in make_indices(dims[:-1]):
            try:
                index = tuple(index)
            except TypeError:
                pass

            # Sort trace
            sx = np.sort(tx[index])
            # Append to list
            intervals[index] = calc_min_interval(sx, alpha)
        # Transpose back before returning
        return np.array(intervals)
    else:
        # Sort univariate node
        sx = np.sort(x)
        return np.array(calc_min_interval(sx, alpha))
    
def hpd_multiple(x):
    x = np.array(x)
    hpd99 = hpd_single(x, 0.01)
    hpd95 = hpd_single(x, 0.05)
    hpd50 = hpd_single(x, 0.50)
    
    return [hpd99[0], hpd95[0], hpd50[0], np.median(x), hpd50[1], hpd95[1], hpd99[1]]


def toYearFraction(datestring):
    from datetime import datetime as dt
    date = dt.strptime(datestring, '%Y-%m-%d')
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def toYearFraction(datestring):
    import time
    from datetime import datetime as dt
    date = dt.strptime(datestring, '%Y-%m-%d')
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction
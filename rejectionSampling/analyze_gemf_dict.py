#!/usr/bin/env python3
import numpy as np
import os
import subprocess
from collections import OrderedDict
import pickle
import math
import pandas as pd
import numpy as np


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


def date_gemf_dicts(rej_sam_path, gemf_dicts_path, test=False):
    '''
    For each combination in the rejection sampling, assign a rejection sampling dict to the index case time    
    '''
    with (open(gemf_dicts_path, "rb")) as openfile:
        gemf_dicts = pickle.load(openfile)
        
    rej_sam_df = pd.read_csv(rej_sam_path, sep='\t')
    gemf_extended_dicts = dict()
    count = 0
    for line in open(rej_sam_path):
        if 'BEAST' in line:
            continue
        l = line.strip().split('\t')
        BEASTgeneration = l[0]
        favitesSim = l[2]
        tIndexCase = math.ceil((float(l[5]) - 2019)*365)
        gemf_d = gemf_dicts[favitesSim]
        S0 = gemf_d[0]['S'] + 1 # number of susceptible individuals before pandemic
        gemf_extended_d = dict() # attach GEMF dict to a specific day and then cover the entire 2019-2020 range
        for d in range(1, 365*2 + 1):
            if d < tIndexCase:
                gemf_extended_d[d] = {'S': S0,
                                  'E': 0,
                                  'P1': 0,
                                  'P2': 0,
                                  'I1': 0,
                                  'I2': 0,
                                  'A1': 0,
                                  'A2': 0,
                                  'H': 0,
                                  'R': 0}
            elif (d - tIndexCase) in gemf_d.keys():
                gemf_extended_d[d] = gemf_d[d - tIndexCase]
            else:
                gemf_extended_d[d] = gemf_d[max(gemf_d.keys())]
        gemf_extended_dicts[BEASTgeneration] = gemf_extended_d

        if test:
            count += 1
            if count >= 5:
                break
        
    return gemf_extended_dicts


def pool_dated_gemf_dicts(dated_gemf_dicts):
    '''
    Pool all the timed GEMF dictionaries; 'inf' is cumulative infections 
    '''
    # pooled_gemf_dict = {day:{'S': [], 'E': [], 'P': [], 'P1': [], 'P2': [], 'I': [], 'I1': [], 'I2': [], 'A': [], 'A1': [], 'A2': [], 'H': [], 'R': [], 'inf': []} for day in dated_gemf_dicts[list(dated_gemf_dicts.keys())[0]].keys()}
    pooled_gemf_dict = {day:{'S': [], 'E': [], 'P1': [], 'P2': [], 'I1': [], 'I2': [], 'A1': [], 'A2': [], 'H': [], 'R': [], 'inf': []} for day in dated_gemf_dicts[list(dated_gemf_dicts.keys())[0]].keys()}
    pop_size = dated_gemf_dicts[list(dated_gemf_dicts.keys())[0]][1]['S']
    for run in dated_gemf_dicts:
        gemf_d = dated_gemf_dicts[run]
        for day in gemf_d:
            for comp in gemf_d[day]:
                pooled_gemf_dict[day][comp].append(gemf_d[day][comp])
            # pooled_gemf_dict[day]['P'].append(gemf_d[day]['P1'] + gemf_d[day]['P2'])
            # pooled_gemf_dict[day]['I'].append(gemf_d[day]['I1'] + gemf_d[day]['I2'])
            # pooled_gemf_dict[day]['A'].append(gemf_d[day]['A1'] + gemf_d[day]['A2'])
            pooled_gemf_dict[day]['inf'].append(pop_size - gemf_d[day]['S'])
    return pooled_gemf_dict


def pooled_gemf_stats(pooled_gemf_dict):
    '''
    Get statistics on the pooled gemf dict for plotting purposes
    '''
    pooled_gemf_stats_dict = {'day': [],
                                'S': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'E': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                # 'P': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'P1': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'P2': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                # 'I': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'I1': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'I2': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                # 'A': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'A1': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'A2': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'H': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'R': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]},
                                'inf': {0.5:[], 1:[], 2.5:[], 25:[], 50:[], 75:[], 97.5:[], 99:[], 99.5:[]} #cumulative infections
                                }
    
    for day in pooled_gemf_dict:
        pooled_gemf_stats_dict['day'].append(day)
        for comp in pooled_gemf_dict[day]:
            stats = hpd_multiple(pooled_gemf_dict[day][comp])
            pooled_gemf_stats_dict[comp][0.5].append(stats[0])
            pooled_gemf_stats_dict[comp][2.5].append(stats[1])
            pooled_gemf_stats_dict[comp][25].append(stats[2])
            pooled_gemf_stats_dict[comp][50].append(stats[3])
            pooled_gemf_stats_dict[comp][75].append(stats[4])
            pooled_gemf_stats_dict[comp][97.5].append(stats[5])
            pooled_gemf_stats_dict[comp][99.5].append(stats[6])

    return pooled_gemf_stats_dict


# main function
if __name__ == "__main__":
    # parse args
    from sys import stdin,stdout; from gzip import open as gopen; import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', '--gemf', required=True, type=str, help="GEMF dict; output of create_gemf_dict.py")
    parser.add_argument('-r', '--rejection', required=True, type=str, help="rejection sampling results")
    parser.add_argument('-c', '--coupled', required=False, action='store_true', help='combining lineages A and B?')
    parser.add_argument('-r2', '--rejection2', required=False, type=str, help="addition rejection sampling results")
    parser.add_argument('-test', '--test', required=False, action='store_true', help="for testing purposes; just do the first five rows of the rejection sampling results")
    parser.add_argument('-o', '--output_string', required=True, type=str, help="output string/prefix")
    args,unknown = parser.parse_known_args()

    if args.coupled == False:
        gemf_extended_dicts = date_gemf_dicts(args.rejection, args.gemf, args.test)
        pooled_gemf_dict = pool_dated_gemf_dicts(gemf_extended_dicts)
        pooled_gemf_stats_dict = pooled_gemf_stats(pooled_gemf_dict)

        prefix = args.output_string
        if prefix.endswith('.'):
            prefix = prefix.strip('.')
        OUTPUT_gemf_extended = prefix + '.timedGEMF.pickle'
        OUTPUT_gemf_stats = prefix + '.timedGEMF.stats.pickle'

        with open(OUTPUT_gemf_extended, 'wb') as handle:
            pickle.dump(pooled_gemf_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(OUTPUT_gemf_stats, 'wb') as handle:
            pickle.dump(pooled_gemf_stats_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    else:
        gemf_extended_dicts_1 = date_gemf_dicts(args.rejection, args.gemf, args.test)
        pooled_gemf_dict_1 = pool_dated_gemf_dicts(gemf_extended_dicts_1)
        pooled_gemf_stats_dict_1 = pooled_gemf_stats(pooled_gemf_dict_1)

        gemf_extended_dicts_2 = date_gemf_dicts(args.rejection2, args.gemf, args.test)
        pooled_gemf_dict_2 = pool_dated_gemf_dicts(gemf_extended_dicts_2)
        pooled_gemf_stats_dict_2 = pooled_gemf_stats(pooled_gemf_dict_2)

        pooled_gemf_dicts_comb = dict()
        for day in pooled_gemf_dict_1:
            pooled_gemf_dicts_comb[day] = dict()
            for compartment in pooled_gemf_dict_1[day]:
                pooled_gemf_dicts_comb[day][compartment] = np.array(pooled_gemf_dict_1[day][compartment]) + np.array(pooled_gemf_dict_2[day][compartment])
        pooled_gemf_stats_dict_comb = pooled_gemf_stats(pooled_gemf_dicts_comb)

        prefix = args.output_string
        if prefix.endswith('.'):
            prefix = prefix.strip('.')
        OUTPUT_gemf_extended_1 = prefix + '.timedGEMF_1.pickle'
        OUTPUT_gemf_stats_1 = prefix + '.timedGEMF_1.stats.pickle'

        OUTPUT_gemf_extended_2 = prefix + '.timedGEMF_2.pickle'
        OUTPUT_gemf_stats_2 = prefix + '.timedGEMF_2.stats.pickle'

        OUTPUT_gemf_extended_comb = prefix + '.timedGEMF_combined.pickle'
        OUTPUT_gemf_stats_comb = prefix + '.timedGEMF_combined.stats.pickle'

        with open(OUTPUT_gemf_extended_1, 'wb') as handle:
            pickle.dump(pooled_gemf_dict_1, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(OUTPUT_gemf_stats_1, 'wb') as handle:
            pickle.dump(pooled_gemf_stats_dict_1, handle, protocol=pickle.HIGHEST_PROTOCOL)

        with open(OUTPUT_gemf_extended_2, 'wb') as handle:
            pickle.dump(pooled_gemf_dict_2, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(OUTPUT_gemf_stats_2, 'wb') as handle:
            pickle.dump(pooled_gemf_stats_dict_2, handle, protocol=pickle.HIGHEST_PROTOCOL)

        with open(OUTPUT_gemf_extended_comb, 'wb') as handle:
            pickle.dump(pooled_gemf_dicts_comb, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(OUTPUT_gemf_stats_comb, 'wb') as handle:
            pickle.dump(pooled_gemf_stats_dict_comb, handle, protocol=pickle.HIGHEST_PROTOCOL)





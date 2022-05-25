#!/usr/bin/env python3
from datetime import datetime
from glob import glob
from os import chdir, getcwd, makedirs
from os.path import isdir, isfile
from random import sample, uniform
from subprocess import call
from sys import argv, stdout
from treeswift import read_tree_newick
import argparse

# useful constants
VERSION = '0.0.1'
LOGFILE = None
C_INT_MAX = 2147483647
GEMF_PARA_FN = 'para.txt'
GEMF_NETWORK_FN = 'network.txt'
GEMF_STATUS_FN = 'status.txt'
GEMF_OUTPUT_FN = 'output.txt'
GEMF_LOG_FN = 'log.txt'
# SEQGEN_LOG_FN = 'seqgen_log.txt'
SAAPPHIIRE_PARAMS = [
    # transition rates
    's_to_e_seed', 'e_to_p1', 'p1_to_p2', 'p2_to_i1', 'p2_to_a1', 'i1_to_i2', 'i1_to_h', 'i1_to_r', 'i2_to_h', 'i2_to_r',
    'a1_to_a2', 'a2_to_r', 'h_to_r', 's_to_e_by_e', 's_to_e_by_p1', 's_to_e_by_p2', 's_to_e_by_i1', 's_to_e_by_i2',
    's_to_e_by_a1', 's_to_e_by_a2',

    # initial state frequencies
    'freq_s', 'freq_e', 'freq_p1', 'freq_p2', 'freq_i1', 'freq_i2', 'freq_a1', 'freq_a2', 'freq_h', 'freq_r',
]

# check if user is just printing version
if '--version' in argv:
    print("FAVITES-COVID-Lite version %s" % VERSION); exit()

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print to the log (None implies stdout only)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    print(tmp, end=end); stdout.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()

# Let L = L0, L1, L2, ... be rates, and let X = X0, X1, X2, ... be exponential random variables where Xi has rate Li. Return the probability that Xi is the minimum of X
def prob_exp_min(i, L):
    assert i >= 0 and i < len(L), "Invalid index i. Must be 0 <= i < |L|"
    return L[i]/sum(L)

# roll a weighted die (keys = faces, values = probabilities)
def roll(orig_die):
    assert len(orig_die) != 0, "Empty weighted die"
    total = float(sum(orig_die.values()))
    die = {face:orig_die[face]/total for face in orig_die}
    faces = sorted(die.keys())
    probs = [die[key] for key in faces]
    cdf = [probs[0]]
    while len(cdf) < len(probs):
        cdf.append(cdf[-1] + probs[len(cdf)])
    num = uniform(0, 1)
    index = 0
    while cdf[index] < num:
        index += 1
    return faces[index]

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
    parser.add_argument('--cn_n', required=True, type=int, help="Contact Network: BA Parameter 'n' (number of nodes)")
    parser.add_argument('--cn_m', required=True, type=int, help="Contact Network: BA Parameter 'm' (number of edges from new node)")
    for tn_param in SAAPPHIIRE_PARAMS:
        parser.add_argument('--tn_%s' % tn_param, required=True, type=float, help="Transmission Network: SAAPPHIIRE Parameter '%s'" % tn_param)
    parser.add_argument('--tn_num_seeds', required=True, type=int, help="Transmission Network: SAAPPHIIRE Parameter: Number of Seeds")
    parser.add_argument('--tn_end_time', required=True, type=float, help="Transmission Network: SAAPPHIIRE Parameter: End Time")
    parser.add_argument('--pt_eff_pop_size', required=True, type=float, help="Phylogeny (Time): Coalescent Intra-host Effective Population Size")
    parser.add_argument('--pm_mut_rate', required=True, type=float, help="Phylogeny (Mutations): Mutation Rate")
    parser.add_argument('--path_ngg_barabasi_albert', required=False, type=str, default='ngg_barabasi_albert', help="Path to 'ngg_barabasi_albert' executable")
    parser.add_argument('--path_gemf', required=False, type=str, default='GEMF', help="Path to 'GEMF' executable")
    parser.add_argument('--path_coatran_constant', required=False, type=str, default='coatran_constant', help="Path to 'coatran_constant' executable")
    parser.add_argument('--path_seqgen', required=False, type=str, default='seq-gen', help="Path to 'seq-gen' executable")
    parser.add_argument('--gzip_output', action='store_true', help="Gzip Compress Output Files")
    parser.add_argument('--version', action='store_true', help="Display Version")
    return parser.parse_args()

# simulate contact network under BA model
def simulate_contact_network_ba(n, m, out_fn, path_ngg_barabasi_albert):
    print_log("Contact Network Model: Barabasi-Albert (BA)")
    for param in ['n', 'm']:
        print_log("BA Parameter '%s': %d" % (param, locals()[param]))
    command = [path_ngg_barabasi_albert, str(n), str(m)]
    print_log("NiemaGraphGen Command: %s" % ' '.join(command))
    out_file = open(out_fn, 'w')
    if call(command, stdout=out_file) != 0:
        raise RuntimeError("Error when running: %s" % path_ngg_barabasi_albert)
    out_file.close()
    print_log("Contact Network simulation complete: %s" % out_fn)

# simulate transmission network under SAAPPHIIRE model
def simulate_transmission_network_saapphiire(
        # transition rates
        s_to_e_seed, e_to_p1, p1_to_p2, p2_to_i1, p2_to_a1, i1_to_i2, i1_to_h, i1_to_r, i2_to_h, i2_to_r,
        a1_to_a2, a2_to_r, h_to_r, s_to_e_by_e, s_to_e_by_p1, s_to_e_by_p2, s_to_e_by_i1, s_to_e_by_i2,
        s_to_e_by_a1, s_to_e_by_a2,

        # initial state frequencies
        freq_s, freq_e, freq_p1, freq_p2, freq_i1, freq_i2, freq_a1, freq_a2, freq_h, freq_r,

        # other simulation parameters
        num_seeds, end_time,

        # paths
        cn_fn, gemf_tmp_dir, out_fn, path_gemf):

    # print initial log info
    print_log("Transmission Network Model: SAAPPHIIRE")
    for param in SAAPPHIIRE_PARAMS:
        print_log("SAAPPHIIRE Parameter '%s': %s" % (param, locals()[param]))
    makedirs(gemf_tmp_dir, exist_ok=True)

    # initialize GEMF stuff
    gemf_state_to_num = {'S':0, 'E':1, 'P1':2, 'P2':3, 'I1':4, 'I2':5, 'A1':6, 'A2':7, 'H':8, 'R':9}
    global gemf_num_to_state; gemf_num_to_state = {gemf_state_to_num[state]:state for state in gemf_state_to_num}
    freq_sum = 0
    for state in gemf_state_to_num.keys():
        param = "freq_%s" % state.lower()
        freq = locals()[param]
        if freq < 0:
            raise ValueError("'%s' must be at least 0" % param)
        freq_sum += freq
    if abs(freq_sum - 1) > 0.000001:
        raise ValueError("Sum of SAAPPHIIRE state frequencies must equal 1")

    # create GEMF para.txt file: https://github.com/niemasd/FAVITES/blob/master/modules/TransmissionTimeSample_SAAPPHIIREGEMF.py#L68-L102
    para_fn = "%s/%s" % (gemf_tmp_dir, GEMF_PARA_FN)
    para_file = open(para_fn, 'w')
    para_file.write("[NODAL_TRAN_MATRIX]\n")
    para_file.write("0\t" + str(s_to_e_seed) + "\t0\t0\t0\t0\t0\t0\t0\t0\n")
    para_file.write("0\t0\t" + str(e_to_p1) + "\t0\t0\t0\t0\t0\t0\t0\n")
    para_file.write("0\t0\t0\t" + str(p1_to_p2) + "\t0\t0\t0\t0\t0\t0\n")
    para_file.write("0\t0\t0\t0\t" + str(p2_to_i1) + "\t0\t" + str(p2_to_a1) + "\t0\t0\t0\n")
    para_file.write("0\t0\t0\t0\t0\t" + str(i1_to_i2) + "\t0\t0\t" + str(i1_to_h) + "\t" + str(i1_to_r) + "\n")
    para_file.write("0\t0\t0\t0\t0\t0\t0\t0\t" + str(i2_to_h) + "\t" + str(i2_to_r) + "\n")
    para_file.write("0\t0\t0\t0\t0\t0\t0\t" + str(a1_to_a2) + "\t0\t0\n")
    para_file.write("0\t0\t0\t0\t0\t0\t0\t0\t0\t" + str(a2_to_r) + "\n")
    para_file.write("0\t0\t0\t0\t0\t0\t0\t0\t0\t" + str(h_to_r) + "\n")
    para_file.write("0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n")
    para_file.write("\n[EDGED_TRAN_MATRIX]\n")
    para_file.write("0\t" + str(s_to_e_by_e)  + "\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n\n")
    para_file.write("0\t" + str(s_to_e_by_p1) + "\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n\n")
    para_file.write("0\t" + str(s_to_e_by_p2) + "\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n\n")
    para_file.write("0\t" + str(s_to_e_by_i1) + "\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n\n")
    para_file.write("0\t" + str(s_to_e_by_i2) + "\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n\n")
    para_file.write("0\t" + str(s_to_e_by_a1) + "\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n\n")
    para_file.write("0\t" + str(s_to_e_by_a2) + "\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n\n")
    para_file.write("[STATUS_BEGIN]\n0\n\n")
    infectious = ['E', 'P1','P2','I1','I2','A1','A2']
    para_file.write("[INDUCER_LIST]\n" + ' '.join([str(gemf_state_to_num[i]) for i in infectious]) + "\n\n")
    para_file.write("[SIM_ROUNDS]\n1\n\n")
    para_file.write("[INTERVAL_NUM]\n1\n\n")
    para_file.write("[MAX_TIME]\n" + str(end_time) + "\n\n")
    para_file.write("[MAX_EVENTS]\n" + str(C_INT_MAX) + "\n\n")
    para_file.write("[DIRECTED]\n0\n\n")
    para_file.write("[SHOW_INDUCER]\n1\n\n")
    para_file.write("[DATA_FILE]\n" + '\n'.join([GEMF_NETWORK_FN]*len(infectious)) + "\n\n")
    para_file.write("[STATUS_FILE]\n%s\n\n" % GEMF_STATUS_FN)
    para_file.write("[OUT_FILE]\n%s" % GEMF_OUTPUT_FN)
    para_file.close()
    print_log("Wrote GEMF '%s' file: %s" % (GEMF_PARA_FN, para_fn))

    # write GEMF network file
    network_fn = "%s/%s" % (gemf_tmp_dir, GEMF_NETWORK_FN)
    network_file = open(network_fn, 'w')
    max_node_label = -1
    for l in open(cn_fn):
        if l.startswith("EDGE"):
            parts = l.split('\t'); u = int(parts[1]); v = int(parts[2])
            network_file.write("%d\t%d\n" % (u,v))
            if u > max_node_label:
                max_node_label = u
            if v > max_node_label:
                max_node_label = v
    network_file.close()
    print_log("Wrote GEMF '%s' file: %s" % (GEMF_NETWORK_FN, network_fn))

    # write GEMF status file
    out_file = open(out_fn, 'w')
    status_fn = "%s/%s" % (gemf_tmp_dir, GEMF_STATUS_FN)
    status_file = open(status_fn, 'w')
    seeds = set(sample(range(max_node_label+1), k=num_seeds))
    gemf_state = list()
    for node in range(max_node_label+1):
        if node in seeds:
            status_file.write("%d\n" % gemf_state_to_num['P1'])
            gemf_state.append(gemf_state_to_num['P1'])
            out_file.write("None\t%d\t0\n" % node)
        else:
            status_file.write("%d\n" % gemf_state_to_num['S'])
            gemf_state.append(gemf_state_to_num['S'])
    status_file.close()
    print_log("Wrote GEMF '%s' file: %s" % (GEMF_STATUS_FN, status_fn))

    # run GEMF
    gemf_log_file = open("%s/%s" % (gemf_tmp_dir, GEMF_LOG_FN), 'w')
    orig_dir = getcwd()
    chdir(gemf_tmp_dir)
    command = [path_gemf]
    print_log("GEMF Command: %s" % ' '.join(command))
    if call(command, stdout=gemf_log_file) != 0:
        raise RuntimeError("Error when running: %s" % path_gemf)
    gemf_log_file.close()
    chdir(orig_dir)
    print_log("Finished running GEMF")

    # reload edge-based matrices for ease of use
    matrices = open("%s/%s" % (gemf_tmp_dir, GEMF_PARA_FN)).read().strip()
    outside_infection_matrix = [[float(e) for e in l.split()] for l in matrices[matrices.index('[NODAL_TRAN_MATRIX]'):matrices.index('\n\n[EDGED_TRAN_MATRIX]')].replace('[NODAL_TRAN_MATRIX]\n','').splitlines()]
    matrices = [[[float(e) for e in l.split()] for l in m.splitlines()] for m in matrices[matrices.index('[EDGED_TRAN_MATRIX]'):matrices.index('\n\n[STATUS_BEGIN]')].replace('[EDGED_TRAN_MATRIX]\n','').split('\n\n')]
    matrices = {gemf_state_to_num[infectious[i]]:matrices[i] for i in range(len(infectious))}
    matrices[gemf_state_to_num['S']] = outside_infection_matrix

    # convert GEMF output to FAVITES transmission network format
    for line in open("%s/%s" % (gemf_tmp_dir, GEMF_OUTPUT_FN)):
        t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.split()]
        v,pre,post = int(v),int(pre),int(post)
        lists = lists.split('],[')
        lists[0] += ']'
        lists[-1] = '[' + lists[-1]
        for i in range(1,len(lists)-1):
            if '[' not in lists[i]:
                lists[i] = '[' + lists[i] + ']'
        lists = [eval(l) for l in lists]
        us = []
        for l in lists:
            us.extend(l)
        if post == gemf_state_to_num['R']:
            out_file.write("%s\t%s\t%s\n" % (v,v,t))
        elif gemf_num_to_state[pre] == 'S' and gemf_num_to_state[post] == 'E':
            uRates = [matrices[gemf_state[u]][pre][post] for u in us]
            die = {us[i]:prob_exp_min(i, uRates) for i in range(len(us))}
            if len(die) != 0:
                u = roll(die) # roll die weighted by exponential infectious rates
            elif len(die) == 0 or u == v: # new seed
                u = "None"
            out_file.write("%s\t%s\t%s\n" % (u,v,t))
        gemf_state[v] = post
    print_log("Transmission Network simulation complete: %s" % out_fn)

# sample people the first time they're ascertained
def sample_first_ascertained(gemf_output_fn, out_fn):
    print_log("Sample Time Model: Fist Time Ascertained")
    sample_time = dict()
    for line in open(gemf_output_fn):
        t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.split()]
        v = int(v); post = int(post)
        if v not in sample_time and gemf_num_to_state[post] in {'I1','I2'}:
            sample_time[v] = float(t)
    out_file = open(out_fn, 'w')
    for node in sample_time:
        out_file.write("%d\t%s\n" % (node, sample_time[node]))
    out_file.close()
    print_log("Sample time selection complete: %s" % out_fn)

# simulate phylogeny (unit of time) under coalescent with constant intra-host effective population size
def simulate_phylogeny_coalescent_constant(eff_pop_size, tn_fn, st_fn, out_fn, path_coatran_constant):
    print_log("Phylogenetic Model: Coalescent with Constant Effective Population Size")
    print_log("Coalescent Parameter 'effective population size': %s" % eff_pop_size)
    command = [path_coatran_constant, tn_fn, st_fn, str(eff_pop_size)]
    print_log("CoaTran Command: %s" % ' '.join(command))
    out_file = open(out_fn, 'w')
    if call(command, stdout=out_file) != 0:
        raise RuntimeError("Error when running: %s" % path_coatran_constant)
    out_file.close()
    print_log("Phylogeny simulation complete: %s" % out_fn)

# scale phylogeny (to unit of mutations) using constant mutation rate
def scale_tree_mutation_rate_constant(mut_rate, pt_fn, out_fn):
    print_log("Mutation Rate Model: Constant")
    print_log("Mutation Rate: %s" % mut_rate)
    time_tree = read_tree_newick(pt_fn)
    if isinstance(time_tree, list):
        raise RuntimeError("No phylogeny was simulated. Try changing your parameters or running again")
    for node in time_tree.traverse_preorder():
        if node.edge_length is not None:
            node.edge_length *= mut_rate
        if not node.is_leaf() and node.label is not None: # seq-gen doesn't support internal node labels
            node.label = None
    time_tree.write_tree_newick(out_fn)
    print_log("Mutation rate scaling complete: %s" % out_fn)

# main execution
if __name__ == "__main__":
    # parse and check user args
    args = parse_args()
    if isdir(args.output) or isfile(args.output):
        raise ValueError("Output directory exists: %s" % args.output)
    if args.cn_n < args.cn_m:
        raise ValueError("In the BA model, the number of nodes (%d) must be larger than the number of edges from new nodes (%d)" % (args.cn_n, args.cn_m))
    if args.tn_num_seeds < 1:
        raise ValueError("Must have at least 1 seed individual")

    # set up output directory
    makedirs(args.output, exist_ok=True)
    LOGFILE_fn = "%s/log.txt" % args.output
    LOGFILE = open(LOGFILE_fn, 'w')

    # print run info to log
    print_log("===== RUN INFORMATION =====")
    print_log("FAVITES-COVID-Lite Version: %s" % VERSION)
    print_log("FAVITES-COVID-Lite Command: %s" % ' '.join(argv))
    print_log("Output Directory: %s" % args.output)
    print_log()

    # simulate contact network
    cn_fn = "%s/contact_network.txt" % args.output
    print_log("===== CONTACT NETWORK =====")
    simulate_contact_network_ba(args.cn_n, args.cn_m, cn_fn, args.path_ngg_barabasi_albert)
    print_log()

    # simulate transmission network
    tn_fn = "%s/transmission_network.txt" % args.output
    tn_gemf_tmp_dir = "%s/GEMF_files" % args.output
    print_log("===== TRANSMISSION NETWORK =====")
    simulate_transmission_network_saapphiire(
        # transition rates
        args.tn_s_to_e_seed, args.tn_e_to_p1, args.tn_p1_to_p2, args.tn_p2_to_i1, args.tn_p2_to_a1, args.tn_i1_to_i2, args.tn_i1_to_h, args.tn_i1_to_r, args.tn_i2_to_h, args.tn_i2_to_r,
        args.tn_a1_to_a2, args.tn_a2_to_r, args.tn_h_to_r, args.tn_s_to_e_by_e, args.tn_s_to_e_by_p1, args.tn_s_to_e_by_p2, args.tn_s_to_e_by_i1, args.tn_s_to_e_by_i2,
        args.tn_s_to_e_by_a1, args.tn_s_to_e_by_a2,

        # initial state frequencies
        args.tn_freq_s, args.tn_freq_e, args.tn_freq_p1, args.tn_freq_p2, args.tn_freq_i1, args.tn_freq_i2, args.tn_freq_a1, args.tn_freq_a2, args.tn_freq_h, args.tn_freq_r,

        # other simulation parameters
        args.tn_num_seeds, args.tn_end_time,

        # paths
        cn_fn, tn_gemf_tmp_dir, tn_fn, args.path_gemf)
    print_log()

    # determine sample times
    st_fn = "%s/sample_times.txt" % args.output
    print_log("===== SAMPLE TIMES =====")
    sample_first_ascertained("%s/%s" % (tn_gemf_tmp_dir, GEMF_OUTPUT_FN), st_fn)
    print_log()
    
    # simulate phylogeny (unit of time)
    pt_fn = "%s/tree.time.nwk" % args.output
    print_log("===== PHYLOGENY =====")
    simulate_phylogeny_coalescent_constant(args.pt_eff_pop_size, tn_fn, st_fn, pt_fn, args.path_coatran_constant)
    print_log()

    # scale phylogeny (to unit of mutations)
    pm_fn = "%s/tree.mutations.nwk" % args.output
    print_log("===== MUTATION RATE =====")
    scale_tree_mutation_rate_constant(args.pm_mut_rate, pt_fn, pm_fn)
    print_log()

    # gzip output files (if requested)
    if args.gzip_output:
        print_log("===== GZIP OUTPUT FILES =====")
        fns_to_compress = [cn_fn, tn_fn, st_fn, pt_fn, pm_fn] + list(glob('%s/*' % tn_gemf_tmp_dir))
        for fn in fns_to_compress:
            command = ['gzip', '-9', fn]
            if call(command) != 0:
                raise RuntimeError("Failed to gzip: %s" % fn)
            print_log("Successfully compressed: %s" % fn)

    # clean things up
    LOGFILE.close()
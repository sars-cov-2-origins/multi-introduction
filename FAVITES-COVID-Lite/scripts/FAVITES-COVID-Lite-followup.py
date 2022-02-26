#!/usr/bin/env python3
from datetime import datetime
from glob import glob
from os import chdir, getcwd, makedirs, listdir
from os.path import isdir, isfile
from random import sample, uniform
from subprocess import call
from sys import argv, stdout
from treeswift import read_tree_newick
import argparse
from gzip import open as gopen 

# command: ../../../scripts/FAVITES-COVID-Lite_followup.py -i . -e 0.16438 --pt_eff_pop_size 1 --pm_mut_rate 0.0008 --path_coatran_constant /usr/local/bin/coatran_constant --gzip_output
# command: /path/to/FAVITES-COVID-Lite_followup.py -i /path/to/inputDir -e 0.16438 --pt_eff_pop_size 1 --pm_mut_rate GET_FROM_NEW_RESULTS --path_coatran_constant /path/to/coatran_constant --gzip_output
# need to write a function to subsample ascertained individuals...maybe 
# parallel --jobs 8 ~/scripts/FAVITES-COVID-Lite_followup.py -i {} -e 0.191781 --pt_eff_pop_size 1 --pm_mut_rate 0.00092005 --path_coatran_constant coatran_constant --gzip_output ::: $(seq -w 1098 1100)

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


# print to the log (None implies stdout only)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    print(tmp, end=end); stdout.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()


def subsample_tn(tn_path, out_fn, num_sample=50000):
    if tn_path.lower().endswith('.gz'):
        tn = gopen(tn_path)
    else:
        tn = open(tn_path)

    count = 0
    out_file = open(out_fn, 'w')
    infected_persons = set()
    for l in tn:
        if isinstance(l,bytes):
            u,v,t = l.decode().strip().split('\t')
        else:
            u,v,t = l.strip().split('\t')

        if u == 'None':
            infected_persons.add(v)
            out_file.write('%s\t%s\t%s\n' % (u,v,t))
            count += 1
        elif u != v and count < num_sample:
            infected_persons.add(v)
            out_file.write('%s\t%s\t%s\n' % (u,v,t))
            count += 1
        elif u == v and u in infected_persons:
            out_file.write('%s\t%s\t%s\n' % (u,v,t))

        # if count >= num_sample:
        #     break
    return infected_persons
    # print(count)

# sample people the first time they're ascertained, up to 50k people
def sample_first_50k(gemf_path, out_fn, end_time, method, tn_infected_persons):
    print_log("Sample Time Model: First 50k")
    gemf_state_to_num = {'S':0, 'E':1, 'P1':2, 'P2':3, 'I1':4, 'I2':5, 'A1':6, 'A2':7, 'H':8, 'R':9}
    gemf_num_to_state = {gemf_state_to_num[state]:state for state in gemf_state_to_num}
    if end_time > 1: # need this to be in years
        end_time = end_time/365
    gemf_time = dict()
    sample_time = dict()
    if gemf_path.endswith('gz'):
        gemf_output = gopen(gemf_path)
    else:
        gemf_output = open(gemf_path)


    for index, line in enumerate(gemf_output):
        if gemf_path.endswith('gz'):
            t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.decode().split()]
        else:
            t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.split()]
        
        if v not in tn_infected_persons: # only care about first X number of infections, based on the subsampled transmission network
            continue

        v = int(v); post = int(post)

        if index == 0: # keep track of this; have to include
            index_case = v
        if v not in gemf_time:
            gemf_time[v] = {gemf_num_to_state[post]: float(t)}
        else:
            gemf_time[v][gemf_num_to_state[post]] = float(t)

    if method == 1: # sample between case status and end, or presymptomatic and end if no case status, or exposed and end if not yet presymptomatic
        for node in gemf_time:
            if 'I1' in gemf_time[node] or 'I2' in gemf_time[node]:
                if 'R' in gemf_time[node]:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['I1'], gemf_time[node]['R'])
                    except:
                        sample_time[node] = uniform(gemf_time[node]['I2'], gemf_time[node]['R'])
                else:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['I1'], end_time)
                    except:
                        sample_time[node] = uniform(gemf_time[node]['I2'], end_time)

            elif 'A1' in gemf_time[node] or 'A2' in gemf_time[node]:
                if 'R' in gemf_time[node]:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['A1'], gemf_time[node]['R'])
                    except:
                        sample_time[node] = uniform(gemf_time[node]['A2'], gemf_time[node]['R'])
                else:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['A1'], end_time)
                    except:
                        sample_time[node] = uniform(gemf_time[node]['A2'], end_time)

            elif 'P1' in gemf_time[node] or 'P2' in gemf_time[node]:
                if 'R' in gemf_time[node]:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['P1'], gemf_time[node]['R'])
                    except:
                        sample_time[node] = uniform(gemf_time[node]['P2'], gemf_time[node]['R'])
                else:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['P1'], end_time)
                    except:
                        sample_time[node] = uniform(gemf_time[node]['P2'], end_time)
            else:
                sample_time[node] = uniform(gemf_time[node]['E'], end_time)


    elif method == 2: # sample index case and then all ascertained cases, with ascertained sampled between case status and end
        for node in gemf_time:
            # print(node, gemf_time[node])
            if node == index_case: # either sample directly at hosp, or between case status and end
                if 'H' in gemf_time[node]:
                    sample_time[node] = gemf_time[node]['H']
                    continue

                if 'R' in gemf_time[node]:
                    end = gemf_time[node]['R']
                else:
                    end = end_time
                
                if 'I2' in gemf_time[node]:
                    start = gemf_time[node]['I2']
                elif 'I1' in gemf_time[node]:
                    start = gemf_time[node]['I1']
                elif 'A2' in gemf_time[node]:
                    start = gemf_time[node]['A2']
                elif 'A1' in gemf_time[node]:
                    start = gemf_time[node]['A1']
                elif 'P2' in gemf_time[node]:
                    start = gemf_time[node]['P2']
                elif 'P1' in gemf_time[node]:
                    start = gemf_time[node]['P1']
                else:
                    start = gemf_time[node]['E']
                sample_time[node] = uniform(start, end)

            elif 'I1' in gemf_time[node] or 'I2' in gemf_time[node]:
                if 'R' in gemf_time[node]:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['I1'], gemf_time[node]['R'])
                    except:
                        sample_time[node] = uniform(gemf_time[node]['I2'], gemf_time[node]['R'])
                else:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['I1'], end_time)
                    except:
                        sample_time[node] = uniform(gemf_time[node]['I2'], end_time)


    elif method == 3: # sample ascertained and unascertained, ascertained at ascertainment, unascertained between presymtomatic and end 
        for node in gemf_time:
            # print(gemf_time[node])
            if 'I1' in gemf_time[node] or 'I2' in gemf_time[node]:
                try:
                    sample_time[node] = gemf_time[node]['I1']
                except:
                    sample_time[node] = gemf_time[node]['I2']
            elif 'P1' in gemf_time[node] or 'P2' in gemf_time[node]:
                if 'R' in gemf_time[node]:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['P1'], gemf_time[node]['R'])
                    except:
                        sample_time[node] = uniform(gemf_time[node]['P2'], gemf_time[node]['R'])
                else:
                    try:
                        sample_time[node] = uniform(gemf_time[node]['P1'], end_time)
                    except:
                        sample_time[node] = uniform(gemf_time[node]['P2'], end_time)
            else:
                sample_time[node] = uniform(gemf_time[node]['E'], end_time)

    elif method == 4:
        '''
        Sample index case, but besides the index case, 
        only sample after we hit the first hospitalized individual. Then, 
        we sample between ascertained and recovery
        '''
        first_hosp = False
        first_hosp_node = None
        first_hosp_time = 999999
        for node in gemf_time:
            if 'H' in gemf_time[node]:
                hosp_time = gemf_time[node]['H']
                # print(gemf_time[node]['H'] - gemf_time[node]['I1'])
                if hosp_time < first_hosp_time:
                    first_hosp_node = node
                    first_hosp_time = hosp_time
                    print(first_hosp_time)

        for node in gemf_time:
            # print(node, gemf_time[node])
            if node == index_case: # either sample directly at hosp, or between case status and end
                if 'I1' in gemf_time[node]:
                    start = gemf_time[node]['I1']
                elif 'A1' in gemf_time[node]:
                    start = gemf_time[node]['A1']
                elif 'P1' in gemf_time[node]:
                    start = gemf_time[node]['P1']
                else:
                    start = gemf_time[node]['E']

                if 'R' in gemf_time[node]:
                    end = gemf_time[node]['R']
                else:
                    end = end_time
            
                # if end_time > first_hosp_time and start < first_hosp_time: # if sampling, can't sample before the earliest hosp.
                #     start = first_hosp_time

                sample_time[node] = {'sample_time' : uniform(start, end), 'hosp' : False}
                if 'H' in gemf_time[node]:
                    if sample_time[node]['sample_time'] >= gemf_time[node]['H']:
                        sample_time[node]['hosp'] = True

            elif 'I1' in gemf_time[node] or 'I2' in gemf_time[node]:
                start = gemf_time[node]['I1']
                if 'R' in gemf_time[node]: # specify end of sampling time window
                    end = gemf_time[node]['R']
                else:
                    end = end_time

                if 'H' in gemf_time and start < first_hosp_time: # all hospitalized individuals sampled, so make sure that they're not accidentally sampled before the first hospitalization time 
                    start = first_hosp_time

                sample_time[node] = {'sample_time' : uniform(start, end), 'hosp' : False}
                if 'H' in gemf_time[node]: # is this a hospital sample? 
                    if sample_time[node]['sample_time'] >= gemf_time[node]['H']:
                        sample_time[node]['hosp'] = True



    out_file = open(out_fn, 'w')
    # count = 0
    print(len(sample_time))
    if method < 4:
        for node in sample_time:
            out_file.write("%d\t%s\n" % (node, sample_time[node]))
            # count += 1
            # if count == 50000:
            #     break
    else:
        '''
        1. sort sample_time by sample_time
        2. Start iterating through the dictionary
        3. Write the sample time of the index case
        4. Write the sample time of the first hospitalized case and then 
           keep writing sample times until we have written 50k sample times 
           (or reached the end of the samples)
        '''
        sample_time_sorted = dict(sorted(sample_time.items(), key=lambda item: item[1]['sample_time']))
        out_file.write("%d\t%s\n" % (index_case, sample_time[index_case]['sample_time']))
        # count = 1
        # print(first_hosp_time)
        for node in sample_time_sorted:
            if node == index_case or sample_time_sorted[node]['sample_time'] < first_hosp_time:
                continue
            out_file.write("%d\t%s\n" % (node, sample_time_sorted[node]['sample_time']))
            # count += 1
            # if count == 1000:
            #     break


    out_file.close()
    print_log("Sample time selection complete: %s" % out_fn)


# simulate phylogeny (unit of time) under coalescent with constant intra-host effective population size
def simulate_phylogeny_coalescent_constant(eff_pop_size, tn_subsampled_fn, st_fn, out_fn, path_coatran_constant):
    print_log("Phylogenetic Model: Coalescent with Constant Effective Population Size")
    print_log("Coalescent Parameter 'effective population size': %s" % eff_pop_size)
    command = [path_coatran_constant, tn_subsampled_fn, st_fn, str(eff_pop_size)]
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


# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-tn', '--transmission_network', required=True, type=str, help="Transmission network")
    parser.add_argument('-g', '--gemf', required=True, type=str, help="GEMF output file")
    parser.add_argument('-e', '--end_time', required=True, type=float, help="Simulation end time")
    parser.add_argument('--pt_eff_pop_size', required=True, type=float, help="Phylogeny (Time): Coalescent Intra-host Effective Population Size")
    parser.add_argument('--pm_mut_rate', required=True, type=float, help="Phylogeny (Mutations): Mutation Rate")
    parser.add_argument('--path_coatran_constant', required=True, type=str, default='coatran_constant', help="Path to 'coatran_constant' executable")
    parser.add_argument('-m', '--method', required=True, type=float, help="Choose {1,2,3,4}: (1) sample all (between A/I and R); (2) sample index (at hosp or between case and R) and I (between I and R); (3) sample all (I at I, A between A and R); ")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output directory")
    parser.add_argument('--gzip_output', action='store_true', help="Gzip Compress Output Files")
    return parser.parse_args()


# main execution
if __name__ == "__main__":
    # parse and check user args
    args = parse_args()
    # if not isdir(args.output) and not isfile(args.output):
    #     raise ValueError("Input directory does not exist: %s" % args.output)

    # for path in listdir(args.output):
    #     if 'transmission' in path:
    #         tn_path = args.output + '/' + path
    #     elif 'GEMF_files' in path:
    #         for gemf_path in listdir(args.output + '/' + path):
    #             if 'output' in gemf_path:
    #                 gemf_output_path = args.output + '/' + path + '/' + gemf_path
    tn_path = args.transmission_network
    gemf_output_path = args.gemf
    # print(tn_path, gemf_output_path)

    # set up output directory
    # makedirs(args.output, exist_ok=True)
    LOGFILE_fn = "%s/log_followup.txt" % args.output
    LOGFILE = open(LOGFILE_fn, 'w')

    # print run info to log
    print_log("===== RUN INFORMATION =====")
    print_log("Output Directory: %s" % args.output)
    print_log()

    # subsample transmission network
    print_log("===== TRANSMISSION SUBSAMPLING =====")
    print_log("Output Directory: %s" % args.output)
    # tn_subsampled_fn = tn_path.split('txt')[0] + 'subsampled.txt'
    tn_subsampled_fn = '%s/transmission_network.subsampled.txt' % args.output
    tn_infected_persons = subsample_tn(tn_path, tn_subsampled_fn, num_sample=50000)
    print_log()
  
    # determine sample times
    st_fn = "%s/subsample_times.txt" % args.output
    print_log("===== SAMPLE TIMES =====")
    sample_first_50k(gemf_output_path, st_fn, args.end_time, args.method, tn_infected_persons)
    print_log()
    
    # simulate phylogeny (unit of time)
    pt_fn = "%s/tree.time.subsampled.nwk" % args.output
    print_log("===== PHYLOGENY =====")
    # print(tn_subsampled_fn, st_fn, pt_fn)
    simulate_phylogeny_coalescent_constant(args.pt_eff_pop_size, tn_subsampled_fn, st_fn, pt_fn, args.path_coatran_constant)
    print_log()

    # scale phylogeny (to unit of mutations)
    pm_fn = "%s/tree.mutations.subsampled.nwk" % args.output
    print_log("===== MUTATION RATE =====")
    scale_tree_mutation_rate_constant(args.pm_mut_rate, pt_fn, pm_fn)
    print_log()

    # gzip output files (if requested)
    if args.gzip_output:
        print_log("===== GZIP OUTPUT FILES =====")
        fns_to_compress = [tn_subsampled_fn, st_fn, pt_fn, pm_fn] 
        print(fns_to_compress)
        for fn in fns_to_compress:
            command = ['gzip', '-9', fn]
            if call(command) != 0:
                raise RuntimeError("Failed to gzip: %s" % fn)
            print_log("Successfully compressed: %s" % fn)

    # clean things up
    LOGFILE.close()

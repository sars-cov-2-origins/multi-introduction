#!/usr/bin/env python3
'''
Pull data from all GEMF and coalescence files from the simulations

Assumed very specific data structure
'''
import treeswift
import os
import subprocess
import string
import math

# main function
if __name__ == "__main__":
    # parse args
    from sys import stdin,stdout; from gzip import open as gopen; import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-dir', '--dir', required=True, type=str, help="simulation directory")
    parser.add_argument('-d', '--days', required=True, type=int, help="days in simulation")
    args,unknown = parser.parse_known_args()

    results_d = {}
    for run in sorted(os.listdir(args.dir)):
        if set(run.lower()).intersection(set(string.ascii_lowercase)):
            continue
        coal_path = args.dir + '/' + run + '/coalData_parameterized.txt'
        gemf_path = args.dir + '/' + run + '/GEMF_files/firsts.txt'

        # stable coalescence
        results_d[run] = {}
        coal_results = [x.strip().split('\t') for x in open(coal_path).readlines()]
        stable_coal = coal_results[-1][1]
        results_d[run]['TimeToCoal'] = float(stable_coal)
        results_d[run]['infections'] = coal_results[(math.ceil(float(stable_coal)*365) + 1)][2]
        # number infections by stable coalescence
        # num_inf = coal_results[-1][2]
        # for coal_day in coal_results[::-1]:
        #     if coal_day[1] != stable_coal:
        #         break
        #     num_inf = coal_day[2]
        # results_d[run]['infections'] = num_inf

        # for line in open(coal_path):
        #     if 'time' in line.lower():
        #         continue
        #     l = line.strip().split('\t')
        #     coal_time = float(l[1])
        #     results_d[run]['TimeToCoal'] = coal_time

        for line in open(gemf_path):
            # print(run, line)
            l = line.strip().split('\t')
            try:
                results_d[run][l[0]] = float(l[1])
            except:
                results_d[run][l[0]] = args.days/365

    OUTPUT = args.dir + '/FAVITES_results.tsv'
    with open(OUTPUT, 'w') as f:
        f.write('Run\tCoalescent time\tFirst ascertained\tFirst unascertained\tFirst hospitalized\tinfections\n')
        for run in results_d:
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (run, results_d[run]['TimeToCoal'], results_d[run]['I'], results_d[run]['A'], results_d[run]['H'], results_d[run]['infections']))



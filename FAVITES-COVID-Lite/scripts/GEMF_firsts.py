#!/usr/bin/env python3
'''
Get first occurrences of ascertained, unascertained, and hospitalized persons

Just putting it in the same directory as GEMF files
'''
import treeswift
import os
import subprocess
import numpy as np
import pandas as pd
import math

# main function
if __name__ == "__main__":
    # parse args
    from sys import stdin,stdout; from gzip import open as gopen; import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', '--gemf', required=True, type=str, help="GEMF output")
    args,unknown = parser.parse_known_args()

    OUTPUT = '/'.join(args.gemf.split('/')[:-1] + ['firsts.txt'])
    gemf_state_to_num = {'S':0, 'E':1, 'P1':2, 'P2':3, 'I1':4, 'I2':5, 'A1':6, 'A2':7, 'H':8, 'R':9}
    gemf_num_to_state = {gemf_state_to_num[state]:state for state in gemf_state_to_num}

    firsts = {'I': None, 'A': None, 'H': None}
    if args.gemf.endswith('gz'):
        gemf_output = gopen(args.gemf)
    else:
        gemf_output = open(args.gemf)

    for index, line in enumerate(gemf_output):
        if args.gemf.endswith('gz'):
            t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.decode().split()]
        else:
            t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.split()]
        post = gemf_num_to_state[int(post)]
        t = float(t)
        
        if (post == 'I1' or post == 'I2') and firsts['I'] == None:
            firsts['I'] = t
        elif (post == 'A1' or post == 'A2') and firsts['A'] == None:
            firsts['A'] = t
        elif post == 'H' and firsts['H'] == None:
            firsts['H'] = t

        if firsts['I'] and firsts['A'] and firsts['H']:
            break
    with open(OUTPUT, 'w') as f:
        for key in firsts:
            f.write('%s\t%s\n' % (key, firsts[key]))



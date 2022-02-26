#!/usr/bin/env python3
from os import chdir, getcwd, makedirs, listdir
from os.path import isdir, isfile
import argparse
import math
from gzip import open as gopen
import pickle


# create dict from GEMF run
def gemf_dict(gemf_path):
    gemf_d = dict()
    gemf_cum_d = dict()
    day = 0
    count = 0

    if gemf_path.endswith('.gz'):
        gemf = gopen(gemf_path)
    else:
        gemf = open(gemf_path)

    for index, line in enumerate(gemf):
        if gemf_path.endswith('gz'):
            t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.decode().split()]
        else:
            t,rate,v,pre,post,num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,lists = [i.strip() for i in line.split()]
        
        day = math.ceil(float(t)*365)
        post = int(post); num0 = int(num0); num1 = int(num1); num2 = int(num2); num3 = int(num3); num4 = int(num4); num5 = int(num5); num6 = int(num6); num7 = int(num7); num8 = int(num8); num9 = int(num9)
        
        if index == 0:
            gemf_d[day] = {'S' : num0 , 
                         'E' : num1 , 
                         'P1' : num2 , 
                         'P2' : num3 ,
                         'I1' : num4 , 
                         'I2' : num5 , 
                         'A1' : num6 ,
                         'A2' : num7 , 
                         'H' : num8 , 
                         'R' : num9 
                         }

            gemf_cum_d[day] = {'S' : num0 , 
                             'E' : 1 , 
                             'P1' : 1 , 
                             'P2' : 1 ,
                             'I1' : num4 , 
                             'I2' : num5 , 
                             'A1' : num6 ,
                             'A2' : num7 , 
                             'H' : num8 , 
                             'R' : num9 
                             }
            old_day = day

        # if it's a new day, make an entry for the new day in the dicts
        elif day not in gemf_d:
            # print(day, post, num0, num1, num2, num3, num4, num5, num6, num7, num8, num9)
            gemf_cum_d[day] = gemf_cum_d[old_day].copy()
            gemf_cum_d[day][gemf_num_to_state[post]] += 1
            gemf_cum_d[day]['S'] = num0

            gemf_d[day] = {'S' : num0 , 
                         'E' : num1 , 
                         'P1' : num2 , 
                         'P2' : num3 ,
                         'I1' : num4 , 
                         'I2' : num5 , 
                         'A1' : num6 ,
                         'A2' : num7 , 
                         'H' : num8 , 
                         'R' : num9 
                            }
            old_day = day

        # if it's the same day, update the dicts
        else:
            gemf_cum_d[day][gemf_num_to_state[post]] += 1
            gemf_cum_d[day]['S'] = num0
            gemf_d[day] = {'S' : num0 , 
                         'E' : num1 , 
                         'P1' : num2 , 
                         'P2' : num3 ,
                         'I1' : num4 , 
                         'I2' : num5 , 
                         'A1' : num6 ,
                         'A2' : num7 , 
                         'H' : num8 , 
                         'R' : num9 
                            }
            old_day = day    

    for day in range(max(gemf_d.keys()) + 1):
        if day not in gemf_d:
            gemf_d[day] = gemf_d[day - 1].copy()
            gemf_cum_d[day] = gemf_cum_d[day - 1].copy()
    return gemf_d, gemf_cum_d


# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('-i', '--input', required=True, type=str, help="gemf file; for testing purposes; delete later")
    parser.add_argument('-d', '--dir', required=True, type=str, help="Directory of simulations; assumes directory/*runs*/GEMF_files/output.txt.gz file structure")
    return parser.parse_args()


 # main execution
if __name__ == "__main__":
    args = parse_args()

    gemf_state_to_num = {'S':0, 'E':1, 'P1':2, 'P2':3, 'I1':4, 'I2':5, 'A1':6, 'A2':7, 'H':8, 'R':9}
    gemf_num_to_state = {gemf_state_to_num[state]:state for state in gemf_state_to_num}

    all_gemf_dicts = {}
    all_gemf_cum_dicts = {}
    for path in sorted(listdir(args.dir)):
        print(path)
        gemf_path = args.dir + '/' + path + '/GEMF_files/output.txt.gz'
        if not isfile(gemf_path):
            gemf_path = gemf_path[:-3]
            if not isfile(gemf_path):
                continue

        gemf_d, gemf_cum_d = gemf_dict(gemf_path)
        gemf_d = {k: gemf_d[k] for k in sorted(gemf_d.keys())}
        gemf_cum_d = {k: gemf_cum_d[k] for k in sorted(gemf_cum_d.keys())}
        all_gemf_dicts[path] = gemf_d
        all_gemf_cum_dicts[path] = gemf_cum_d


    with open('FAVITES_GEMF_dict.pickle', 'wb') as handle:
        pickle.dump(all_gemf_dicts, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('FAVITES_GEMF_cumulative_dict.pickle', 'wb') as handle:
        pickle.dump(all_gemf_cum_dicts, handle, protocol=pickle.HIGHEST_PROTOCOL)



    


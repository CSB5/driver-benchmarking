#!/usr/bin/env python                                                                                             

import sys
sys.path.append('/mnt/software/unstowable/hotnet2')

import argparse
import os
import os.path
import json

parser = argparse.ArgumentParser()
parser.add_argument('delta_dir', help="directory containing delta subdirectories")
parser.add_argument('sig_threshold', help="significance threshold", type=float)
args = parser.parse_args()

# Select smallest delta (minimum edge weight) for each run with the most
# significant (P<sig_threshold) subnetwork sizes k
   
# get all subdirectory names in delta directory starting with 'delta_'
d = args.delta_dir
subdir = [os.path.join(d,o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]
subdir = [i for i in subdir if i.startswith(os.path.join(d,'delta_'))]


# The method I used to chose \delta is based on email correspondance with Max Leiserson [mdml@cs.brown.edu]:
# 1. Choose a (k, P-value) per \delta. We chose the k with the smallest P-value, or returned null if the smallest 
# P-value was > 0.05 (i.e. there was no significance).
# 2. Choose \delta. We ran the above procedure for each \delta to get a (k, P-value), and chose the \delta 
# with the smallest P-value. In the case of ties we took the \delta with the smaller k.


# for all subdirectories, get p-values from significance.txt
# and count the number of p-values <= sig_threshold (n)
deltas = []
k = []
pval = []

for num_dir in range(0,len(subdir)):
    with open((subdir[num_dir]+'/significance.txt'), 'r') as f:
        deltas.append(float(subdir[num_dir].split('delta_')[1]))
        pval_temp = 1
        k_temp = 0
        for l in f.readlines()[1:]:
            l = l.split()
            if (float(l[-1]) < pval_temp): # if p values are same lower k is kept
                pval_temp = float(l[-1])
                k_temp = int(l[0])
        pval.append(pval_temp)
        k.append(k_temp)
 
# chose the \delta with the smallest P-value.                                                                         
minpval_index = [i for i, x in enumerate(pval) if x == min(pval)]
k_min = min([k[index] for index in minpval_index])
# In the case of ties we took the \delta with the smaller k.  
selected_delta_index = [i for i in minpval_index if k[i] == k_min]
delta_max = max([deltas[index] for index in selected_delta_index])

# in case there are two (k, pval) that are the same I select the maximum /delta
# I inferred that this is what the HotNet2 team did by looking at their data and /deltas
# for example: http://compbio-research.cs.brown.edu/pancancer/hotnet2/public/data/runs/by-type/gbm/iref/freq/
if(len(selected_delta_index) > 1): 
    selected_delta_index = [i for i in selected_delta_index if deltas[i] == delta_max]

# option 1: load results.json file for the chosen delta, edit to remove 
# networks that fall below the p-value
import json

with open(subdir[selected_delta_index[0]]+'/results.json') as json_file:
    json_data = json.load(json_file)

for key in json_data['sizes'].keys():
    temp = int(key)
    if (temp > 10): temp = 10
    if (json_data['statistics'][str(temp)]['pval'] >= args.sig_threshold):
        num_networks = len(json_data['components'])
        for j in range(1,num_networks+1):
            if (temp <10):
                if (len(json_data['components'][(num_networks-j)]) == temp):
                    del json_data['components'][(num_networks-j)]
            else:
                if (len(json_data['components'][(num_networks-j)]) >= 10):
                    del json_data['components'][(num_networks-j)]  
        del json_data['sizes'][key] 
        
with open(os.path.join(args.delta_dir, 'results.json'), 'w') as outfile:
    json.dump(json_data, outfile, sort_keys=True, indent=4, separators=(',', ': '))
# option 2: copy results.json file for the chosen delta as is to delta_dir
#import shutil
#shutil.copy2(os.path.join(subdir[nmax_index], 'results.json'), 
#             os.path.join(args.delta_dir, 'results.json'))



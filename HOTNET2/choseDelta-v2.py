#!/usr/bin/env python                                                                                             

import sys
sys.path.append('/mnt/software/unstowable/hotnet2')

import argparse
import os
import os.path
import json
import numpy as np
from itertools import izip

parser = argparse.ArgumentParser()
parser.add_argument('delta_dir', help="directory containing delta subdirectories")
parser.add_argument('freq_file', help="mutation frequency data file")
parser.add_argument('output_dir', help="directory to save .result file in")
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

if(pval[selected_delta_index[0]] <= args.sig_threshold):
    # find maximum network size which still contains significant networks
    significance_table = np.loadtxt(subdir[selected_delta_index[0]]+'/significance.txt', skiprows = 1)
    K = k[selected_delta_index[0]]-3
    while K >= 0:
        if(significance_table[K,2]-significance_table[K+1,2] > 0):
            if(significance_table[K,3] > args.sig_threshold):
                break
            else:
                K = K-1
        else:
            K = K-1

    genes = []

    # get all genes in networks with size > significance_table[K,0]
    with open(subdir[selected_delta_index[0]]+'/components.txt','r') as f:
        for l in f.readlines():
            l = l.split()
            if len(l) <= significance_table[K,0]: break
            genes.extend(l)
        
    genes.sort()

    # ranked genes based on mutation frequency
    scores = [0]*len(genes)
    i = 0
    with open(args.freq_file,'r') as f:
        for l in f.readlines():
            l = l.split()
            if l[0] in genes: 
                scores[genes.index(l[0])] = l[1]
         
    data = sorted(izip(genes, scores), reverse=True, key=lambda x: x[1])
    rank = range(1,len(scores)+1)


    # write to .result file
    if not os.path.exists(args.output_dir): os.makedirs(args.output_dir)

    f = open(os.path.join(args.output_dir,'HotNet2.result'), 'w')
    f.write('Gene_name\tSample\tRank\tScore\tInfo\n')
    for i in range(0,len(genes)):
        cola, colc = data[i]
        f.write(cola + '\tALL\t' + str(rank[i]) + '\t' + str(colc) + '\n')
    f.close()
else: # no deltas/networks greater than significance threshold
    f = open(os.path.join(args.output_dir,'HotNet2.result'), 'w')
    f.write('Gene_name\tSample\tRank\tScore\tInfo\n')
    f.close()


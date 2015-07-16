#!/usr/bin/env python  

import os
import os.path
import re
from itertools import izip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('consensus_file', help="file containing consensus networks")
parser.add_argument('heat_file', help="list of genes and corresponding heat scores")
parser.add_argument('output_dir', help="directory to store output files")
args = parser.parse_args()

genes = []
l_genes = []         

# get line number to start reading in genes
for i, l in enumerate(open(args.consensus_file, 'r')):
    if re.search('#Consensus', l): 
        l_genes = (i+1)
        break
    
# get all gene names in core consensus network (in square brackets)
with open(args.consensus_file,'r') as f:
    for l in f.readlines()[l_genes:]:
        text = re.search('\[(.*?)\]', l)
        found = []
        if text:
            found = text.group(1) 
            found = found.replace(' ','')
            found = found.split(',')
        genes.extend(found)
        
genes.sort()

# read in heat scores corresponding to the gene list
# Note: genes in heat score file must be in alphabetical order!
scores = []
i = 0
with open(args.heat_file,'r') as f:
    for l in f.readlines():
        if (i < len(genes)):
            if re.search(genes[i], l):
                i = i+1
                l = l.split()
                scores.append(l[1])

# rank the genes based on the heat scores
data = sorted(izip(genes, scores), reverse=True, key=lambda x: x[1])
rank = range(1,len(scores)+1)

# write to .result file
if not os.path.exists(args.output_dir): os.makedirs(args.output_dir)
        
f = open(os.path.join(args.output_dir,'HotNet2.result'), 'w')
f.write('Gene_name\tSample\tRank\tScore\tInfo\n')
for i in range(0,len(data)):
    cola, colc = data[i]
    f.write(cola + '\t\t' + str(rank[i]) + '\t' + colc + '\n')
f.close()


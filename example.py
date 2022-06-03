import os
from collections import Counter
import sys
sys.setrecursionlimit(3000)
import pandas as pd
import debruijn as db

sequences = []
sequences_scores = []
for root, dir, files in os.walk('avastin/avastin'):
    root = root + '/'
    for file in files:
        filename = root + file
        data = pd.read_csv(filename, delimiter='\t')
        temp = data[data['Score']>0.8]
        temp = temp[-50<temp['PPM Difference']]
        temp = temp[temp['PPM Difference']<50]
        temp.reset_index(inplace=True)
        sequences.extend(temp['DENOVO'].values)
        for i in range(len(temp)):
            sequences_scores.append([temp['DENOVO'][i],temp['Score'][i]])
            sequences.append(temp['DENOVO'][i])


sequences = Counter(sequences)
print(len(sequences))
print(sequences)

sequences = list(sequences.keys())

k = 3
g = db.construct_graph(sequences, k)
# print_graph(g)
# for k in g.keys():
#   print k, g[k]
# g = construct_graph(reads)
contig = db.output_contigs(g)
print(contig)

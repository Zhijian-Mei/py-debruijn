import os
from collections import Counter
import sys
import pandas as pd
import debruijn_new as db

sequences = []
sequences_scores = []
for root, dir, files in os.walk('avastin/avastin'):
    root = root + '/'
    for file in files:
        filename = root + file
        data = pd.read_csv(filename, delimiter='\t')
        temp = data[data['Score']>0.6]
        temp = temp[-50<temp['PPM Difference']]
        temp = temp[temp['PPM Difference']<50]
        temp.reset_index(inplace=True)
        sequences.extend(temp['DENOVO'].values)
        for i in range(len(temp)):
            sequences_scores.append([temp['DENOVO'][i],temp['Score'][i]])
            sequences.append(temp['DENOVO'][i])


sequences = Counter(sequences)
sequences = list(sequences.keys())
print(len(sequences))


k = 4
step = 50
contig = []
for i in range(0,len(sequences),step):
    try:
        input = sequences[i:i+step]
    except:
        input = sequences[i:]
    print('number of kmer in this step: ',len(input))
    g = db.construct_graph(input, k)
    output = db.output_contigs(g)
    if len(contig) == 0:
        contig.extend(output)
    else:
        for j in contig:
            for m in output:
                if j[-k:] == m[:k]:
                    j = j + m[k:]
                    output.remove(m)
                elif j[:k] == m[-k:]:
                    j = m + j[k:]
                    output.remove(m)
        contig.extend(output)
for item in contig:
    print(item,len(item))

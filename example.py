import json
import os

from collections import Counter
import sys

sys.setrecursionlimit(10000)
import pandas as pd
import test_debruijn as db


def getScore(edge_count_table, contig, k):
    score = 0
    for i in range(len(contig) - k):
        score += edge_count_table[contig[i:i + k + 1]]
    return score


sequences = []
score_cut = 0.5
threshold = 1
k_lowerlimit = 5
k_upperlimit = 10

for root, dir, files in os.walk('Ab_1/Ab_1'):
    root = root + '/'
    for file in files:
        filename = root + file
        data = pd.read_csv(filename, delimiter='\t')
        temp = data[data['Score'] >= score_cut]
        temp = temp[-50<temp['PPM Difference']]
        temp = temp[temp['PPM Difference']<50]
        sequences.extend(temp['DENOVO'].values)
        temp.reset_index(inplace=True)
        sequences.extend(temp['DENOVO'].values)
        for i in range(len(temp)):
            sequences.append(temp['DENOVO'][i])

sequences = Counter(sequences)
sequences = list(sequences.keys())
print(len(sequences))


for k in range(k_lowerlimit, k_upperlimit + 1):
    if k <= k_upperlimit - 1:
        g, pull_out_read, branch_kmer, already_pull_out, edge_count_table = db.construct_graph(sequences, k,
                                                                                               threshold=threshold)
    else:
        g, pull_out_read, branch_kmer, already_pull_out, edge_count_table = db.construct_graph(sequences, k,
                                                                                               threshold=threshold, final=True)

    sequences = db.output_contigs(g, branch_kmer, already_pull_out)
    sequences.sort(key=lambda x: getScore(edge_count_table, x, k), reverse=True)
    if k == k_upperlimit:
        froot = 'Ab_1_{}-{}mer_{}_{}'.format(k_lowerlimit,k_upperlimit, score_cut,threshold)
        os.mkdir(froot)
        setting = {'score_cut': score_cut, 'threshold': threshold, 'k_lowerlimit': k_lowerlimit,
                   'k_upperlimit': k_upperlimit}
        with open(f'{froot}/setting.json','w') as fw:
            json.dump(setting,fw,indent=4)
        outFile = open(f'{froot}/{froot}.fasta', mode='a+')
        for i in range(len(sequences)):
            outFile.writelines('>SEQUENCE_{}_{}mer\n{}\n'.format(i, k, sequences[i]))
        outFile.close()
    print('max length: ', len(max(sequences, key=lambda x: len(x))))
    print('number of output for k={}: '.format(k), len(sequences))
    if k <= k_upperlimit - 1:
        sequences.extend(pull_out_read)
        print('number of pull out read: ', len(pull_out_read))


from pprint import pprint

import pandas as pd
from tqdm import trange


def read_fasta(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    dic = {}

    for i in range(len(lines)):
        line = lines[i]
        if line[0] == '>':
            id = line.split(' ')[0][1:]
            contig = lines[i + 1]
            dic[id] = contig.rstrip()
    return dic


def findOverlap(intervals):
    v1 = max([pair[0] for pair in intervals])
    v2 = min([pair[1] for pair in intervals])
    if v1 > v2:
        return None
    return [v1, v2]


template_name = 'templates/homo_template.fasta'
froot = 'avastin_5-10mer_0.6_2'
df = pd.read_csv(f'{froot}/{froot}_blasthomoTemplate.m8', delimiter='\t', header=None)
template_dic = read_fasta(template_name)

df = df[df[2] >= 80]
df = df.sort_values(by=1)
df = df.reset_index(drop=True)

previous_protein = None
protein_sequences_records = {}
for i in trange(len(df)):
    current_protein = df[1][i]
    current_contig_id = df[0][i]
    right = df[9][i]
    left = df[8][i]
    if current_protein not in protein_sequences_records.keys():
        protein_sequences_records[current_protein] = {}
    protein_sequences_records[current_protein][current_contig_id] = [left, right]

keys = list(protein_sequences_records.keys())
for i in trange(len(keys)):
    key = keys[i]
    protein_sequences_records[key] = dict(sorted(protein_sequences_records[key].items(), key=lambda x: x[1]))
    sequence_interval_records = protein_sequences_records[key]
    intervals = []
    for sub_key in sequence_interval_records.keys():
        intervals.append(sequence_interval_records[sub_key])
    print(intervals)
    quit()

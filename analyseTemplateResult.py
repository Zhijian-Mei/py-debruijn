import copy

import pandas as pd
from tqdm import trange


def findOverlap(intervals):
    if len(intervals) <= 1:
        return None
    v1 = max([pair[0] for pair in intervals])
    v2 = min([pair[1] for pair in intervals])
    if v1 > v2:
        return None
    return [v1,v2]


froot = 'avastin_5-10mer_0.6_2'
df = pd.read_csv(f'{froot}/{froot}_blasthomoTemplate.m8', delimiter='\t', header=None)

df = df[df[3] >= 80]

df = df.reset_index(drop=True)

coverage_record = {}

previous = None
for i in trange(len(df)):
    sequence_label = df[0][i]
    protein_label = df[1][i]
    if sequence_label not in coverage_record.keys():  # new sequence
        protein_record = {}
        coverage_record[sequence_label] = 0
    right = df[9][i]
    left = df[8][i]
    coverage = right - left + 1
    coverage_record[sequence_label] += coverage
    if protein_label not in protein_record.keys():
        protein_record[protein_label] = [[left, right]]
    else:
        before = copy.deepcopy(protein_record[protein_label])
        after = copy.deepcopy(protein_record[protein_label])
        after.append([left, right])

        protein_record[protein_label] = after

        overlap_before = findOverlap(before)
        overlap_after = findOverlap(after)
        print()
        print(before)
        print(after)
        print(overlap_before)
        print(overlap_after)
        if overlap_after != overlap_before:
            coverage_record[sequence_label] -= (overlap_after[1] - overlap_after[0])

print(coverage_record)
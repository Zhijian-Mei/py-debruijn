import copy

import pandas as pd
from tqdm import trange


def findOverlap(intervals):
    N = len(intervals)
    if N <= 1:
        return None
    # First interval
    l = intervals[0][0]
    r = intervals[0][1]

    # Check rest of the intervals
    # and find the intersection
    for i in range(1, N):

        # If no intersection exists
        if intervals[i][0] > r or intervals[i][1] < l:
            return None

        # Else update the intersection
        else:
            l = max(l, intervals[i][0])
            r = min(r, intervals[i][1])
    return [l,r]



froot = 'avastin_5-10mer_0.6_2'
df = pd.read_csv(f'{froot}/{froot}_blasthomoTemplate.m8', delimiter='\t', header=None)

df = df[df[3] >= 80]

df = df.reset_index(drop=True)

coverage_record = {}

previous = None
for i in trange(len(df)):
    sequence_label = df[0][i]
    protein_label = df[1][i]
    if sequence_label not in coverage_record.keys(): # new sequence
        protein_record = {}
        coverage_record[sequence_label] = 0
    right = df[9][i]
    left = df[8][i]
    coverage = right - left + 1
    coverage_record[sequence_label] += coverage
    if protein_label not in protein_record.keys():
        protein_record[protein_label] = [[left,right]]
    else:
        before = copy.deepcopy(protein_record[protein_label])
        after = copy.deepcopy(protein_record[protein_label])
        after.append([left,right])
        overlap_before = findOverlap(before)
        overlap_after = findOverlap(after)
        if overlap_after != overlap_before:
            print(before)
            print(after)
            print(overlap_before)
            print(overlap_after)
            coverage_record[sequence_label] -= (overlap_after[1] - overlap_after[0])

            quit()




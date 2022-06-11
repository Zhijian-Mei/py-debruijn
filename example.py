import os
from collections import Counter
import sys
sys.setrecursionlimit(3000)
import pandas as pd
import test_debruijn as db

sequences = []
sequences_scores = []
for root, dir, files in os.walk('BSA/all'):
    root = root + '/'
    for file in files:
        filename = root + file
        print(filename)
        data = pd.read_csv(filename, delimiter='\t')
        temp = data[data['Score']>=0.1]
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

# sequences = ['EVQLVESGGGLVQPGGSLRLSCAASGYTFTNYGMNWVRQAPGKGLEWVGWLNTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYPHYYGSSHWYFDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYLCNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMLSRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPLEKTLSKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDLAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK']

k_lowerlimit = 5
k_upperlimit = 6
for k in range(k_lowerlimit,k_upperlimit+1):
    print('number of input for k={}'.format(k), len(sequences))
    g, pull_out_read = db.construct_graph(sequences, k, threshold=2)
    sequences = db.output_contigs(g)
    sequences.sort(key=lambda x: len(x))
    # sequences = [sequence for sequence in sequences if len(sequence) > 200]
    print('number of output has length > 200: ',len(sequences))
    for item in sequences:
        if 'MKWV' in item:
            print(item,len(item))
    if k <= k_upperlimit-1:
        sequences.extend(pull_out_read)
        print('number of pull out read: ',len(pull_out_read))



contig = sequences
contig.sort(key=lambda x:len(x))
# for item in contig:
    # if 'MKWVT' in item:
    # print(item,len(item))


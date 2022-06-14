import os

from collections import Counter
import sys
sys.setrecursionlimit(10000)
import pandas as pd
import test_debruijn as db

sequences = []
sequences_scores = []
for root, dir, files in os.walk('BSA/all'):
    root = root + '/'
    for file in files:
        filename = root + file
        data = pd.read_csv(filename, delimiter='\t')
        # temp = data
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
k_upperlimit = 10
for k in range(k_lowerlimit,k_upperlimit+1):
    if k <= k_upperlimit-1:
        g, pull_out_read,pull_out_kmer = db.construct_graph(sequences, k, threshold=2)
    else:
        g, pull_out_read, pull_out_kmer = db.construct_graph(sequences, k, threshold=2, final=True)
    sequences = db.output_contigs(g,pull_out_kmer)
    sequences.sort(key=lambda x: len(x))
    outFile = open('BSA_{}mer.fasta'.format(k),mode='a+')
    for i in range(len(sequences)):
        outFile.writelines('>SEQUENCE_{}\n{}\n'.format(i,sequences[i]))
    outFile.close()
    print('max length: ',len(sequences[-1]))
    print('number of output for k={}: '.format(k),len(sequences))
    if k <= k_upperlimit-1:
        sequences.extend(pull_out_read)
        print('number of pull out read: ',len(pull_out_read))





# print('*' * 100)
# for item in contig:
#     print(item,len(item))

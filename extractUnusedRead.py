import copy
import json
import os
from collections import Counter

import pandas as pd
from generateTemplatesBlastReport import read_fasta
from Bio.Blast.Applications import NcbiblastpCommandline

if __name__ == '__main__':
    froot = 'avastin_5-10mer_0.6_2'
    contig_filepath = f'{froot}/{froot}_modified_sorted.fasta'

    settingFile = open(f'{froot}/setting.json', 'r')
    setting = json.load(settingFile)
    souceFilePath = 'avastin/avastin'
    DF = pd.DataFrame()
    for root, dir, files in os.walk(souceFilePath):
        root = root + '/'
        for file in files:
            filename = root + file
            data = pd.read_csv(filename, delimiter='\t')
            temp = data[data['Score'] >= setting['score_cut']]
            temp = temp[-50 < temp['PPM Difference']]
            temp = temp[temp['PPM Difference'] < 50]
            temp.reset_index(inplace=True, drop=True)
            DF = DF.append(temp)
    DF.reset_index(inplace=True, drop=True)


    contig_dic = read_fasta(contig_filepath)
    contigs = list(contig_dic.values())

    reads = DF['DENOVO'].values
    reads_copy = copy.deepcopy(reads)
    titles = DF['TITLE'].values
    unused_reads = []
    for i in range(len(reads)):
        read = reads[i]
        # for j in range(len(contigs)):
        #     contig = contigs[j]
        #     if read in contig:
        #         break
        #     if read not in contig and (j == len(contigs)-1):
        #         unused_reads.append(read)
        unused_reads.append(read)
    print(Counter(unused_reads))
    quit()
    unused_reads = list(Counter(unused_reads))
    print(unused_reads)
    quit()

    with open(f'{froot}/unusedReads.fasta','w') as f:
        for i in range(len(unused_reads)):
            f.write('>unused_reads_{}\n'.format(i))
            f.write('{}\n'.format(unused_reads[i]))




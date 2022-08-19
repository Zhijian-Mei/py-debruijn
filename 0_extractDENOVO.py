import argparse
import json
import os
from collections import Counter

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-source', type=str, required=True)
    parser.add_argument('-score', type=float, required=True)
    parser.add_argument('-t', type=int, required=True)
    parser.add_argument('-kl', type=int, required=True)
    parser.add_argument('-ku', type=int, required=True)

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    score_cut = args.score
    threshold = args.t
    k_lowerlimit = args.kl
    k_upperlimit = args.ku
    source = f'{args.source}/{args.source}'
    froot = f'{args.source}_{k_lowerlimit}-{k_upperlimit}mer_{score_cut}_{threshold}'

    for root, dir, files in os.walk(source):
        root = root + '/'
        print(root)
        for file in files:
            filename = root + file
            data = pd.read_csv(filename, delimiter='\t')
            temp = data[data['Score'] < score_cut]
            # temp = temp[temp['Score'] > 0]
            unused_reads=temp['DENOVO'].values
            unused_reads = list(Counter(unused_reads).keys())
            for read in unused_reads:
                if type(read) is not str:
                    print(read)
            quit()
            unused_reads = [x for x in unused_reads if len(x) > k_lowerlimit]
            print(unused_reads)
            quit()

    quit()

    df = pd.DataFrame()

    df['Peptide'] = unused_reads
    df['Charge'] = [2 for i in range(len(unused_reads))]
    df['Type'] = ['HCD' for i in range(len(unused_reads))]
    df['NCE'] = [25 for i in range(len(unused_reads))]
    setting = {'score_cut': score_cut, 'threshold': threshold, 'k_lowerlimit': k_lowerlimit,
               'k_upperlimit': k_upperlimit}

    try:
        os.mkdir(froot)
    except:
        pass

    with open(f'{froot}/setting.json', 'w') as fw:
        json.dump(setting, fw, indent=4)
    df.to_csv(f'{froot}/unused_reads.tsv', sep='\t')
    os.system(
        f'python PredFull/predfull.py --input {froot}/unused_reads.tsv --model PredFull/pm.h5 --output {froot}/unused_reads_prediction.mgf'
    )
    os.system(
        f'./msSLASH/bin/bruteforce  -e 1111466_E.mgf -l unused_reads_prediction.mgf -d empty.mgf'
    )

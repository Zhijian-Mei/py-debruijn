import argparse
import os

import pandas as pd

from test_debruijn import read_reads
from generateTemplatesBlastReport import read_fasta

def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    froot = args.froot
    contigs = read_fasta(f'{froot}/{froot}_sorted.fasta')
    template_contig_key = list(contigs.keys())[0]
    template_contig = contigs[template_contig_key]
    contigs.pop(template_contig_key)
    while True:
        template_contig_filename = f'{froot}/{froot}_templateContig.fasta'
        template_contig_filename_temp = f'{froot}/{froot}_templateContig-db.fasta'
        query_filename = f'{froot}/{froot}_queryContigs.fasta'
        with open(template_contig_filename,'w') as f:
            f.write('>template_Contig\n')
            f.write(template_contig)
        with open(query_filename,'w') as f:
            for contigs_key in contigs.keys():
                f.write('>{}\n'.format(contigs_key))
                f.write(contigs[contigs_key])
                f.write('\n')
        os.system(f'prerapsearch -d {template_contig_filename} -n {template_contig_filename_temp}')
        os.system(f'rapsearch -q {query_filename} -d {template_contig_filename_temp} -o {froot}/rapsearch_templateContig -z 6')
        os.system(
            f'python processRapsearchM8.py -input {froot}/rapsearch_templateContig.m8 -output {froot}/rapsearch_templateContig_refactor.m8')
        df = pd.read_csv(f'{froot}/rapsearch_templateContig_refactor.m8', delimiter='\t', header=None)
        print(df)
        quit()

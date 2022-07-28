from Bio.Align.Applications import MafftCommandline
from io import StringIO

from tqdm import trange


def read_fasta(fname):
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    fasta={}

    for i in range(len(lines)):
        line = lines[i]
        if '>' in line:
            fasta[line.rstrip()] = lines[i+1].rstrip()
    return fasta


froot = 'avastin_5-10mer_0.6_2'
contigs = read_fasta(f'{froot}/{froot}_cluster.fasta')
templates = read_fasta(f'templates/homo_templates_final.fasta')
templates_keys = list(templates.keys())
in_file = "input.fasta"
for contig in contigs.keys():
    for i in trange(len(templates_keys)):
        template = templates_keys[i]
        string = contig + '\n' + contigs[contig] + '\n' + template+'\n' + templates[template]
        f = open(in_file, 'w')
        f.write(string)
        f.close()
        mafft_cline = MafftCommandline(input=in_file,clustalout=True)
        stdout, stderr = mafft_cline()
        with open("output.out", "a") as handle:
            handle.write(stdout)
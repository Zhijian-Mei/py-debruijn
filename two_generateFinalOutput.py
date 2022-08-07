import json
import random
from pprint import pprint

from tqdm import trange

from test_debruijn import read_reads


def checkSubSequence(contig, output):
    for o in output:
        if contig in o:
            return False
    return True






if __name__ == '__main__':
    froot = 'avastin_4-10mer_0.6_2'
    contigs = read_reads(f'{froot}/{froot}.fasta')
    with open(f'{froot}/setting.json') as f:
        setting = json.load(f)
    print(setting)
    print('max length before concat: ', len(max(contigs, key=lambda x: len(x))))

    print(len(contigs))
    k = 4
    already_chosen = []
    for i in trange(100000):
        concat = False
        current = random.choice(contigs)
        if current in already_chosen:
            continue
        else:
            already_chosen.append(current)
        for i in range(len(contigs)):
            current_head_to_contig_tail_concat = False
            current_tail_to_contig_head_concat = False
            contig = contigs[i]
            current_head = current[:k]
            contig_tail = contig[len(contig) - k:]
            current_tail = current[len(current)-k:]
            contig_head = contig[:k]
            if current_head == contig_tail: # contig的尾与current的头k个相等
                contigs[i] = contig + current[k:]
                concat = True
                continue
            if current_tail == contig_head: # contig的头与current的尾k个相等
                contigs[i] = current + contig[k:]
                concat = True
                continue

            for j in range(k-1):
                shuffled_current_head = list(current_head)
                shuffled_current_head[j],shuffled_current_head[j+1] = shuffled_current_head[j+1],shuffled_current_head[j]
                shuffled_current_head = ''.join(shuffled_current_head)
                if shuffled_current_head == contig_tail:
                    result_1 = contig + current[k:]
                    result_2 = contig[:len(contig)-k] + current
                    result = random.choice([result_1,result_2])
                    current_head_to_contig_tail_concat = True
                    break
            if current_head_to_contig_tail_concat:
                contigs[i] = result
                concat = True
                continue

            for j in range(k-1):
                shuffled_current_tail = list(current_tail)
                shuffled_current_tail[j],shuffled_current_tail[j+1] = shuffled_current_tail[j+1],shuffled_current_tail[j]
                shuffled_current_tail = ''.join(shuffled_current_tail)
                if shuffled_current_tail == contig_head:
                    result_1 = current + contig[k:]
                    result_2 = current[:len(current)-k] + contig
                    result = random.choice([result_1, result_2])
                    current_tail_to_contig_head_concat = True
                    break

            if current_tail_to_contig_head_concat:
                contigs[i] = result
                concat = True
                continue

        if concat:
            contigs.remove(current)
    print(len(contigs))

    outputs = []
    for contig in contigs:
        if contig not in outputs and checkSubSequence(contig, outputs):
            outputs.append(contig)
    print(len(outputs))
    print('max length after concat: ', len(max(outputs, key=lambda x: len(x))))

    k = setting['k_upperlimit']
    outFile = open(f'{froot}/{froot}_concatenated.fasta', mode='a+')
    for i in range(len(outputs)):
        outFile.writelines('>SEQUENCE_{}_{}mer\n{}\n'.format(i, k, outputs[i]))
    outFile.close()

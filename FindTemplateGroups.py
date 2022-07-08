from pprint import pprint

import pandas as pd

from generateTemplatesBlastReport import read_fasta

if __name__ == '__main__':
    template_name = 'templates/homo_template.fasta'
    froot = 'avastin_5-10mer_0.6_2'
    contig_filepath = f'{froot}/{froot}_modified_sorted.fasta'
    df = pd.read_csv(f'{froot}/{froot}_blasthomoTemplate.m8', delimiter='\t', header=None)
    df = df[df[2] >= 80]
    df = df.sort_values(by=0)
    df = df.reset_index(drop=True)
    template_dic = read_fasta(template_name, 'Homo')
    templates = list(template_dic.keys())
    contig_dic = read_fasta(contig_filepath)
    contigs = list(contig_dic.keys())

    dfList = df.values
    sequence_template_id_pair_dic = {}
    for item in dfList:
        template_id = item[:2][1]
        label = item[:2][0] + '+' + template_id
        value_list = list(item[2:])
        sequence_template_id_pair_dic[label] = value_list

    template_contig_group = {}
    while len(contigs) != 0:
        current_contig = contigs[0]
        template_length_record = {}
        for template_id in templates:
            label = current_contig + '+' + template_id
            try:
                identity = sequence_template_id_pair_dic[label][0]
            except:
                continue
            if identity > 90:
                template_length_record[template_id] = len(template_dic[template_id])
        if len(template_length_record) == 0:
            contigs.remove(current_contig)
            continue
        candidate_templates = sorted(list(template_length_record.items()), key=lambda x: x[1], reverse=True)

        best_template = candidate_templates[0][0]
        # candidates_group = [best_template[0]]
        # for i in range(1,len(candidate_templates)):
        #     if candidate_templates[i][1] == best_template[1]:
        #         candidates_group.append(candidate_templates[i][0])
        #
        # candidate_length_record = {}
        # for candidate in candidates_group:
        #     candidate_length_record[candidate] = len(template_dic[candidate])
        #
        # best_template = sorted(list(candidate_length_record.items()),key=lambda x:x[1],reverse=True)[0][0]


        template_contig_group[best_template] = [current_contig]
        contigs.remove(current_contig)

        remove = []
        for contig in contigs:
            label = contig + '+' + best_template
            try:
                value = sequence_template_id_pair_dic[label]
            except:
                continue
            if value:
                template_contig_group[best_template] += [contig]
                remove.append(contig)
        for item in remove:
            contigs.remove(item)

    pprint(template_contig_group)
    print(template_contig_group.keys())
    print(len(template_contig_group.keys()))

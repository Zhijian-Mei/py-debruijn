import re
from pprint import pprint

import pandas as pd

from generateTemplatesBlastReport import read_fasta

class Template:
    def __init__(self,template_id,template_sequence):
        self.sequence = template_sequence
        self.id = template_id
        self.contigArrays = []
        self.different_position = []

class Contig:
    def __init__(self,contig_id,contig_sequence,template_interval,contig_interval):
        self.sequence = contig_sequence
        self.id = contig_id
        self.template_interval = template_interval
        self.contig_interval = contig_interval

class fillingTemplate:
    def __init__(self,template_sequence):
        self.template_sequence = template_sequence
        self.fill = list(' ' * len(self.template_sequence))


    def fill_match(self,contig):
        contig_sequence = list(contig.sequence[contig.contig_interval[0]-1:contig.contig_interval[1]])
        self.fill[contig.template_interval[0]-1:contig.template_interval[1]] = contig_sequence



    def get_match_result(self):
        return ''.join(self.fill)

def checkOverlap(contig_array, contig):
    intervals = []
    for item in contig_array:
        intervals.append(item.template_interval)
    intervals.append(contig.template_interval)
    intervals.sort(key=lambda x: x[0])

    for i in range(len(intervals) - 1):
        if intervals[i][1] > intervals[i + 1][0]:
            return True
    return False

def toHTML(string):
    result = []
    for char in string:
        if char == ' ':
            result.append('&nbsp;&nbsp;')
        else:
            result.append(char)
    result = ''.join(result)

    return result

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

    # pprint(template_contig_group)

    report_path = f'{froot}/{froot}_TemplateMatchReport.txt'
    outFile = open(report_path, 'w')
    message = ''


    for template_id in template_contig_group.keys():
        template = Template(template_id,template_dic[template_id])
        for contig_id in template_contig_group[template_id]:
            label = contig_id + '+' + template_id
            value = sequence_template_id_pair_dic[label]
            contig = Contig(contig_id,contig_dic[contig_id],[value[6],value[7]],[value[4],value[5]])
            if len(template.contigArrays) > 0:
                overlap = True
                for contig_array in template.contigArrays:
                    if not checkOverlap(contig_array,contig):
                        contig_array.append(contig)
                        overlap = False
                        break
                if overlap:
                    template.contigArrays.append([contig])
            else:
                template.contigArrays.append([contig])


        result_sequences = []
        for contig_array in template.contigArrays:
            contig_array = sorted(contig_array,key=lambda x:x.template_interval[0])
            match_result = fillingTemplate(template.sequence)
            for i in range(len(contig_array)):
                contig = contig_array[i]
                if (contig.contig_interval[1] - contig.contig_interval[0]) != (contig.template_interval[1] - contig.template_interval[0]):
                    continue
                match_result.fill_match(contig)
            result_sequence = match_result.get_match_result()
            if not re.search('[a-zA-Z]', result_sequence):
                continue
            result_sequences.append(result_sequence)
            different_position = [index for index in range(len(template.sequence)) if template.sequence[index] != result_sequence[index] and result_sequence[index] != ' ']
            for position in different_position:
                if position not in template.different_position:
                    template.different_position.append(position)
        message += '-' * 500
        message += '\n'
        message += 'Template ID: {}\n'.format(template.id)
        changed_line = list(' ' * len(template.sequence))
        for position in template.different_position:
            changed_line[position] = '*'
        changed_line = ''.join(changed_line)
        message += changed_line
        message += '\n'
        message += '{}\n'.format(template.sequence)
        for result_sequence in result_sequences:
            message += result_sequence
            message += '\n'


    outFile.write(message)
    outFile.close()


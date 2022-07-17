import json
import os
import re
from pprint import pprint
import ast
import numpy as np


import pandas as pd
from tqdm import trange

from generateTemplatesBlastReport import read_fasta

class Template:
    def __init__(self,template_id,template_sequence):
        self.sequence = template_sequence
        self.id = template_id
        self.contigArrays = []
        self.different_position = []
        self.letters_correctRate = {}
        for i in range(len(self.sequence)):
            self.letters_correctRate[i] = {}

class Contig:
    def __init__(self,contig_id,contig_sequence,template_interval,contig_interval):
        self.sequence = contig_sequence
        self.id = contig_id
        self.template_interval = template_interval
        self.contig_interval = contig_interval
        self.rates = {}
        for i in range(len(self.sequence)):
            self.rates[i] = 1

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

        reads = DF['DENOVO'].values
        position_scores = DF['Positional Score'].values

        for array_index in range(len(template.contigArrays)):
            contig_array = template.contigArrays[array_index]
            contig_array = sorted(contig_array, key=lambda x: x.template_interval[0])
            for contig in contig_array:
                for i in range(len(reads)):
                    read = reads[i]
                    if read in contig.sequence:
                        match = re.search(read,contig.sequence)
                        read_positional_scores = ast.literal_eval(position_scores[i])
                        for j in range(match.start(),match.end()):
                            contig.rates[j] = contig.rates[j] * (1-read_positional_scores[j - match.start()])

        print()
        print('*' * 500)
        print(template.sequence)
        message += '*' * 500
        message += '\n'
        message += 'Template ID: {}'.format(template.id)
        message += '\n'
        message += template.sequence
        message += '\n'
        for contig_array in template.contigArrays:
            contig_array = sorted(contig_array, key=lambda x: x.template_interval[0])
            for contig in contig_array:
                if (contig.contig_interval[1] - contig.contig_interval[0]) != (contig.template_interval[1] - contig.template_interval[0]):
                    continue
                template_points = [x for x in range(contig.template_interval[0]-1,contig.template_interval[1])]
                contig_points = [x for x in range(contig.contig_interval[0] - 1,contig.contig_interval[1])]
                for i in range(len(template_points)):
                    template_point = template_points[i]
                    contig_point = contig_points[i]
                    current_template_position = template.letters_correctRate[template_point]
                    current_contig_letter = contig.sequence[contig_point]
                    if current_contig_letter not in current_template_position.keys():
                        current_template_position[current_contig_letter] = contig.rates[contig_point]
                    else:
                        current_template_position[current_contig_letter] = current_template_position[current_contig_letter] * contig.rates[contig_point]

        position_keys = list(template.letters_correctRate.keys())
        result_sequences = []
        for key in position_keys:
            candidate_letters = template.letters_correctRate[key]
            candidate_letters = dict(sorted(candidate_letters.items(), key=lambda item: item[1]))
            if candidate_letters != {}:
                while len(result_sequences) < len(candidate_letters):
                    result_sequences.append(list(' '*len(template.sequence)))
                letters = list(candidate_letters.keys())
                for i in range(len(letters)):
                    result_sequences[i][key] = letters[i]
        for sequence in result_sequences:
            print(''.join(sequence))
            message += ''.join(sequence)
            message += '\n'



        # result_sequences = []
        # for contig_array in template.contigArrays:
        #     contig_array = sorted(contig_array,key=lambda x:x.template_interval[0])
        #     match_result = fillingTemplate(template.sequence)
        #     for i in range(len(contig_array)):
        #         contig = contig_array[i]
        #         if (contig.contig_interval[1] - contig.contig_interval[0]) != (contig.template_interval[1] - contig.template_interval[0]):
        #             continue
        #         match_result.fill_match(contig)
        #     result_sequence = match_result.get_match_result()
        #     if not re.search('[a-zA-Z]', result_sequence):
        #         continue
        #     result_sequences.append(result_sequence)
        #     different_position = [index for index in range(len(template.sequence)) if template.sequence[index] != result_sequence[index] and result_sequence[index] != ' ']
        #     for position in different_position:
        #         if position not in template.different_position:
        #             template.different_position.append(position)
        # message += '-' * 500
        # message += '\n'
        # message += 'Template ID: {}\n'.format(template.id)
        # changed_line = list(' ' * len(template.sequence))
        # for position in template.different_position:
        #     changed_line[position] = '*'
        # changed_line = ''.join(changed_line)
        # message += changed_line
        # message += '\n'
        # message += '{}\n'.format(template.sequence)
        # for result_sequence in result_sequences:
        #     message += result_sequence
        #     message += '\n'


    outFile.write(message)
    outFile.close()


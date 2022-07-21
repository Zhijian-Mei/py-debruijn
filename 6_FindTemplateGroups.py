import json
import os
import re
from pprint import pprint
import ast
import numpy as np

import pandas as pd
from tqdm import trange
from Bio.Blast.Applications import NcbiblastpCommandline
from generateTemplatesBlastReport import read_fasta
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class Template:
    def __init__(self, template_id, template_sequence):
        self.sequence = template_sequence
        self.id = template_id
        self.contigArrays = []
        self.different_position = []
        self.letters_correctRate = {}
        for i in range(len(self.sequence)):
            self.letters_correctRate[i] = {}
        self.unusedReads_match = {}
        for i in range(len(self.sequence)):
            self.unusedReads_match[i] = []


class Contig:
    def __init__(self, contig_id, contig_sequence, template_interval, contig_interval):
        self.sequence = contig_sequence
        self.id = contig_id
        self.template_interval = template_interval
        self.contig_interval = contig_interval
        self.rates = {}
        for i in range(len(self.sequence)):
            self.rates[i] = 1


class fillingTemplate:
    def __init__(self, template_sequence):
        self.template_sequence = template_sequence
        self.fill = list(' ' * len(self.template_sequence))

    def fill_match(self, contig):
        contig_sequence = list(contig.sequence[contig.contig_interval[0] - 1:contig.contig_interval[1]])
        self.fill[contig.template_interval[0] - 1:contig.template_interval[1]] = contig_sequence

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
            result.append('-')
        else:
            result.append(char)
    result = ''.join(result)

    return result


if __name__ == '__main__':
    template_name = 'templates/homo_templates.fasta'
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

    template_dic = read_fasta(template_name)
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
        best_identity = 0
        best_template = None
        for template_id in templates:
            label = current_contig + '+' + template_id
            try:
                identity = sequence_template_id_pair_dic[label][0]
            except:
                continue
            # if identity > 90:
            #     template_length_record[template_id] = len(template_dic[template_id])
            if identity > best_identity:
                best_identity = identity
                best_template = template_id
        # if len(template_length_record) == 0:
        if not best_template:
            contigs.remove(current_contig)
            continue
        # candidate_templates = sorted(list(template_length_record.items()), key=lambda x: x[1], reverse=True)
        # best_template = candidate_templates[0][0]
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
    html_path = f'{froot}/{froot}_TemplateMatchReport.html'
    htmlFile = open(html_path, 'w')
    html = '''<!DOCTYPE html>
    <body>
    '''
    reads = DF['DENOVO'].values
    position_scores = DF['Positional Score'].values

    for template_id in template_contig_group.keys():
        template = Template(template_id, template_dic[template_id])
        for contig_id in template_contig_group[template_id]:
            label = contig_id + '+' + template_id
            value = sequence_template_id_pair_dic[label]
            contig = Contig(contig_id, contig_dic[contig_id], [value[6], value[7]], [value[4], value[5]])
            if len(template.contigArrays) > 0:
                overlap = True
                for contig_array in template.contigArrays:
                    if not checkOverlap(contig_array, contig):
                        contig_array.append(contig)
                        overlap = False
                        break
                if overlap:
                    template.contigArrays.append([contig])
            else:
                template.contigArrays.append([contig])

        for array_index in range(len(template.contigArrays)):
            contig_array = template.contigArrays[array_index]
            contig_array = sorted(contig_array, key=lambda x: x.template_interval[0])
            for contig in contig_array:
                for i in range(len(reads)):
                    read = reads[i]
                    if read in contig.sequence:
                        match = re.search(read, contig.sequence)
                        read_positional_scores = ast.literal_eval(position_scores[i])
                        for j in range(match.start(), match.end()):
                            contig.rates[j] = contig.rates[j] * (1 - read_positional_scores[j - match.start()])

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
                if (contig.contig_interval[1] - contig.contig_interval[0]) != (
                        contig.template_interval[1] - contig.template_interval[0]):
                    continue
                template_points = [x for x in range(contig.template_interval[0] - 1, contig.template_interval[1])]
                contig_points = [x for x in range(contig.contig_interval[0] - 1, contig.contig_interval[1])]
                for i in range(len(template_points)):
                    template_point = template_points[i]
                    contig_point = contig_points[i]
                    current_template_position = template.letters_correctRate[template_point]
                    current_contig_letter = contig.sequence[contig_point]
                    if current_contig_letter not in current_template_position.keys():
                        current_template_position[current_contig_letter] = contig.rates[contig_point]
                    else:
                        current_template_position[current_contig_letter] = current_template_position[
                                                                               current_contig_letter] * contig.rates[
                                                                               contig_point]

        position_keys = list(template.letters_correctRate.keys())
        result_sequences = []
        for key in position_keys:
            candidate_letters = template.letters_correctRate[key]
            candidate_letters = dict(sorted(candidate_letters.items(), key=lambda item: item[1]))
            if candidate_letters != {}:
                while len(result_sequences) < len(candidate_letters):
                    result_sequences.append(list(' ' * len(template.sequence)))
                letters = list(candidate_letters.keys())
                for i in range(len(letters)):
                    result_sequences[i][key] = letters[i]
        for sequence in result_sequences:
            print(''.join(sequence))
            message += ''.join(sequence)
            message += '\n'

        message += '-' * 100 + '\n'
        print('-' * 100)
        message += 'Unused reads blast result: \n'
        print('Unused reads blast result: ')
        message += template.sequence + '\n'
        print(template.sequence)
        with open(f'{froot}/temp.fasta', 'w') as f:
            f.write('>{}\n'.format(template.id))
            f.write(template.sequence)
        out = f'{froot}/{froot}_unusedReadsBlastTemplate.m8'
        query = f'{froot}/unusedReads.fasta'
        command = NcbiblastpCommandline(query=query,
                                        subject=f'{froot}/temp.fasta',
                                        outfmt=6,
                                        out=out,
                                        )
        command()

        unusedReadsTemplateResults = pd.read_csv(out, delimiter='\t', header=None)
        unusedReadsTemplateResults = unusedReadsTemplateResults[unusedReadsTemplateResults[2] >= 90]
        unusedReadsTemplateResults.reset_index(drop=True, inplace=True)

        unusedReads_dic = read_fasta(query)

        unusedReads_value_dic = {}
        for value in unusedReadsTemplateResults.values:
            unusedReads_value_dic[value[0]] = list(value[1:])

        for unusedRead in unusedReads_value_dic.keys():
            item = unusedReads_value_dic[unusedRead]
            read_left = item[5]
            read_right = item[6]
            template_left = item[7]
            template_right = item[8]
            if (read_right - read_left) != (template_right - template_left):
                continue
            matchedReadSeq = unusedReads_dic[unusedRead][read_left - 1:read_right]
            # if matchedReadSeq == 'SGL':
            #     print(unusedReads_dic[unusedRead])
            #     print(unusedRead)
            #     quit()
            for i in range(template_left - 1, template_right):
                current_read_letter = matchedReadSeq[i - (template_left - 1)]
                if current_read_letter not in template.unusedReads_match[i]:
                    template.unusedReads_match[i] += [current_read_letter]
        unusedReadsResultSequence = []
        for key in template.unusedReads_match.keys():
            while len(template.unusedReads_match[key]) > len(unusedReadsResultSequence):
                unusedReadsResultSequence.append(list(' ' * len(template.sequence)))
            letters = template.unusedReads_match[key]
            for i in range(len(letters)):
                unusedReadsResultSequence[i][key] = letters[i]
        for sequence in unusedReadsResultSequence:
            print(''.join(sequence))
            message += ''.join(sequence) + '\n'

        merged_result = []
        while len(result_sequences) < max(len(result_sequences), len(unusedReadsResultSequence)):
            result_sequences.append(list(' ' * len(template.sequence)))
        for unusedReadResult in unusedReadsResultSequence:
            for i in range(len(template.sequence)):
                letter = unusedReadResult[i]
                for contig_result in result_sequences:
                    if contig_result[i] == letter:
                        break
                    if contig_result[i] == ' ':
                        contig_result[i] = '<font color="green">{}</font>'.format(letter)
                        break

        for sequence in result_sequences:
            for i in range(len(sequence)):
                if sequence[i] != ' ' and len(sequence[i]) == 1:
                    sequence[i] = '<font color="blue">{}</font>'.format(sequence[i])

        merged_result = result_sequences

        step = 250
        print('-' * 100 + 'Merged Result' + '-' * 100)
        html += '*' * 100 + 'Merged Result' + '*' * 100 + '<br>'
        html += 'Template ID: {}<br>'.format(template.id)
        for i in range(0, len(template.sequence), step):
            try:
                sub_template = template.sequence[i:i + step]
            except:
                sub_template = template.sequence[i:]
            print(sub_template)
            html += '<pre>' + sub_template + '</pre>'
            for sequence in merged_result:
                try:
                    sub_sequence = sequence[i:i + step]
                except:
                    sub_sequence = sequence[i:]
                sub_sequence = ''.join(sub_sequence)
                if not any(c.isalpha() for c in sub_sequence):
                    continue
                html += '<pre>' + sub_sequence + '</pre>'

            html += '<br>'

    htmlFile.write(html)
    htmlFile.close()
    outFile.write(message)
    outFile.close()

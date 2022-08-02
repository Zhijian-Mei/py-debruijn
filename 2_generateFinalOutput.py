import json

from test_debruijn import read_reads

def checkSubSequence(contig,output):
    for o in output:
        if contig in o:
            return False
    return True

def checkHeadToTail(c1,contigs,k):
    result = c1
    for c2 in contigs:
        if result[:k] == c2[len(c2)-k:]:
            result = c2 + result[k:]
        elif c2[:k] == result[len(c1) - k:]:
            result = result + c2[k:]
    return result

froot = 'Ab_1_5-10mer_0.6_2'
contigs = read_reads(f'{froot}/{froot}.fasta')
f = open(f'{froot}/setting.json')
setting = json.load(f)
print(setting)
print('max length before concat: ',len(max(contigs,key=lambda x:len(x))))

temp = []
for c1 in contigs:
    result = checkHeadToTail(c1,contigs,3)
    temp.append(result)



output=[]
for contig in temp:
    if contig not in output and checkSubSequence(contig,output):
        output.append(contig)


print('max length after concat: ',len(max(output,key=lambda x:len(x))))

k = setting['k_upperlimit']
outFile = open(f'{froot}/{froot}_modified.fasta', mode='a+')
for i in range(len(output)):
    outFile.writelines('>SEQUENCE_{}_{}mer\n{}\n'.format(i,k, output[i]))
outFile.close()
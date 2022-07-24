
import pandas as pd
from tqdm import trange

froot= 'avastin_5-10mer_0.6_2'
file = open(f'{froot}/rapsearch_outputs.m8',mode='r+')
lines = file.readlines()
file.close()
df = pd.DataFrame(columns=[0,1,2,3,4,5,6,7,8,9,10,11])
outlines = []
aline = []
for i in trange(len(lines)):
    line = lines[i]
    if '#' not in line:
        items = line.split('\t')
        for item in items:
            if item:
                aline.append(item.rstrip())
    if len(aline) == 12:
        outlines.append(aline)
        aline = []

df = pd.DataFrame(outlines,columns=[0,1,2,3,4,5,6,7,8,9,10,11])

df.to_csv(f'{froot}/rapsearch_outputs_refactor.m8',index=False,sep='\t',header=None)
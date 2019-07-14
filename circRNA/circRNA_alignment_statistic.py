import os
import re
import sys
import glob
import pandas as pd

#inputpath = '/data/yaoyh/Anding_hospital/batch_1/star_circ'
inputpath = sys.argv[1]
outfile = sys.argv[2]


mapping_dict = {}
samples = []

for i in next(os.walk(inputpath))[1]:
    samples.append(i)
    for file in glob.glob(os.path.join(inputpath, i, '*final.out')):
        with open(file, 'r') as filereader:
            line_count = 0
            for row in filereader:
                line_count += 1
                row = row.strip()
                item = row.split('|')
                if line_count == 6:
                    mapping_dict.setdefault('Total Input Reads', []).append(item[1].strip())
                elif line_count == 9:
                    mapping_dict.setdefault('Uniquely mapped reads', []).append(item[1].strip())
                elif line_count == 10:
                    mapping_dict.setdefault('Uniquely mapped rate', []).append(item[1].strip())
                elif line_count == 12:
                    mapping_dict.setdefault('Total spliced reads', []).append(item[1].strip())
                elif line_count == 33:
                    mapping_dict.setdefault('circular spliced reads', []).append(item[1].strip())
                else :
                    #print(line_count)
                    continue
print (mapping_dict)
df = pd.DataFrame(data = mapping_dict)
raw_index = {}
for i in range(len(samples)):
    raw_index.setdefault(df.index[i], samples[i])

df = df.rename(index = raw_index)
df.to_csv(outfile)

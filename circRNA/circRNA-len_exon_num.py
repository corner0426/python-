import os
import sys
import glob
#import pybedtools
#from pybedtools import BedTool
import pandas as pd

CIRCexplorer2_path = sys.argv[1]
len_file = sys.argv[2]

f = open(len_file, 'w+')

circ_id = []
for i in next(os.walk(CIRCexplorer2_path))[1]:
	for file in glob.glob(os.path.join(CIRCexplorer2_path, i, '*known.txt')):
		with open(file, 'r') as filereader:
			for line in filereader:
				l = line.strip('\n').split()
				if len(l) < 10: continue
				pos_id = '%s_%s_%s'%(l[0], l[1], l[2])
				if int(l[12]) >= 2 and (pos_id not in circ_id):
					circ_id.append(pos_id)
					if len(l[10].split(',')) == 1:
						f.write('%s\t%s\t%s\t%d\n' % (pos_id, l[9], l[10], int(l[10])))
					else :
						f.write('%s\t%s\t%s\t%d\n' % (pos_id, l[9], l[10], sum([int(i) for i in l[10].split(',')])))
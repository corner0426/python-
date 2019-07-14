#!/usr/bin/env python
#coding: utf-8

import os
import glob
import sys


CIRCexplorer2_path = sys.argv[1]
bed_file = sys.argv[2]
#6列，chr start end name/circ_id score/0 strand

f = open(bed_file, 'w+')

circ_list = []
for i in next(os.walk(CIRCexplorer2_path))[1]:
    for file in glob.glob(os.path.join(CIRCexplorer2_path, i, '*known.txt')):
        with open(file, 'r') as filereader:
            #header = filereader.readline()
            for line in filereader:
                l = line.strip('\n').split()
                if len(l) < 10: continue
                #if int(l[12]) < 2: continue
                circ_id = '%s_%s_%s'%(l[0], l[1], l[2])
                if circ_id in circ_list:
                    continue
                else:
                    circ_list.append(circ_id)
                    f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (l[0], l[1], l[2], circ_id, '0', l[5]))

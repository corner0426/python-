#!/usr/bin/env python
#coding: utf-8

import pybedtools
from pybedtools import BedTool
import sys

circ_bed_file = sys.argv[1]
#circ_bed_file = 'CIRI_6item.bed'
diff_circ_file = sys.argv[2]
#circ_RPM_filter_matrix.csv
ref_bed_file = sys.argv[3]
#ref_bed_file = 'GCF_CriGri_1.0.bed'
circ_anno_file = sys.argv[4]
#circ_anno_file = 'CIRC_anno.txt'
f = open(circ_anno_file, 'w+')

circ_list = []
#f1 = open(diff_circ_file, 'r')
with open(diff_circ_file, 'r') as filereader:
	header = filereader.readline()
	for line in filereader:
		line = line.split(',')
		circ_id = line[0].strip('""')
		circ_list.append(circ_id)
print ('%d length of circ_list' % (len(circ_list)))

ref_bed = pybedtools.BedTool(ref_bed_file)
with open (circ_bed_file, 'r') as filereader:
    for line in filereader:
        line = line.strip('\n')
        l = line.split('\t')
        circ_id = '%s_%s_%s' % (l[0], l[1], l[2])
        if circ_id in circ_list:
            circ_bed = pybedtools.BedTool(line, from_string=True)
            if ref_bed.intersect(circ_bed):
                b_and_a = ref_bed.intersect(circ_bed, F=0.8)
                feature = []
                name = []
                for i in b_and_a:
                    if i[3] not in feature:
                        feature.append(i[3])
                        #feature = feature + ' ' + i[3]
                    if i[4] not in name:
                        name.append(i[4])
                #name = name + ' ' + i[4]

                f.write(line + '\t' + ' '.join(feature) + '\t' + ' '.join(name) + '\n')
            else:
                f.write(line + '\t' + 'NA' + '\t' + 'NA' + '\n')
        else:
            continue

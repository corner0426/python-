#!/usr/bin/env python
#coding: utf-8
#提取表达矩阵，并筛选至少在一半样本中有2个以上back spliced reads

import os
import sys
import glob
import csv

#

#读取shell中制作的sample name文件，按标号数字排序
#/data1/yaoyh/SXJC_circ/CIRCexplorer2_output_files
sample_file = sys.argv[1]

#a diretory contain diretories for each sample
matrix_file = sys.argv[2]

smk_file = open('/data1/yaoyh/SXJC_circ/result/4.Expression/group_only/smk_id.txt', 'r')
smk_id = []
for line in smk_file:
	line = line.strip('\n')
	smk_id.append(line)

non_smk_file = open('/data1/yaoyh/SXJC_circ/result/4.Expression/group_only/non_smk_id.txt', 'r')
non_smk_id = []
for line in non_smk_file:
	line = line.strip('\n')
	non_smk_id.append(line)


files = []
sample_lst = []
if (os.path.isfile(sample_file)):
	try:
		fin = open(sample_file, 'r')
		for line in fin:
			lineLst = line.strip()##是样本名称（包括路径）
			files.append(lineLst)
			#sample_id = lineLst.split('/')[1].split['@'][1]
			#sample_lst.append(sample_id)
	except :
		print ('Error: Can not open sample file')
		exit(1)
else :
	print ("Error: sample information file '%s' not found!" % (sample_file))
	sys.exit(1)
print ('You have %d samples input' %len(files))




circ_dict = {}
circ_id_lst = []

for s in files:
	sample_id = s.split('/')[1].split('@')[1]
	sample_lst.append(sample_id)
	filereader = open(s)
	#header = filereader.readline()
	for line in filereader:
		v = line.split('\t')
		if int(v[12]) < 2: continue  #只保留read数目大于2的环状RNA
		circ_id = v[0] + '_' + v[1] + '_' + v[2]
		if circ_id not in circ_id_lst:
			circ_id_lst.append(circ_id)
		circ_dict.setdefault(circ_id, {})
		circ_dict[circ_id].setdefault(sample_id,int(v[12]))
print ('You have %d circRNAs detected' % (len(circ_id_lst)))

for circ_id in circ_id_lst:
    for x in sample_lst:
        if not x in circ_dict[circ_id]:
            circ_dict[circ_id].setdefault(x, 0)

#circ_remain = 0
with open(matrix_file, 'w+') as csvfile:
    my_writer = csv.DictWriter(csvfile, fieldnames = ["circ_id"] + [x for x in sample_lst])
    my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
    for i in circ_dict:
        circ_dict[i]["circ_id"] = i
        #count = 0
        #for fn in my_writer.fieldnames[1:]:
        #    if circ_dict[i][fn] < int(2) :
        #        count+=1
        #if count <= float(0.5)*len(sample_lst):
        my_writer.writerow(circ_dict[i])
        #    circ_remain += 1
        #else:
        #    continue

#print ('%d circRNA filtered due to low expression among samples' %(len(circ_id_lst)-circ_remain))

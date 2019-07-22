#!/usr/bin/env python
#coding: utf-8
#提取CIRCexploer2基本注释信息，只保留count>2

import os
import sys
import glob
import csv

#

#读取shell中制作的sample name文件，按标号数字排序
#/data1/yaoyh/SXJC_circ/CIRCexplorer2_output_files
sample_file = sys.argv[1]

#a diretory contain diretories for each sample
anno_file = sys.argv[2]

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

#with open(anno_file, 'w+') as f:
f = open(anno_file, 'w+')
f.write('circ_id\thost_gene\texonCount\texonSize\tisoformName\tcircType\n')
for s in files:
	filereader = open(s)
	for line in filereader:
		v = line.split('\t')
		if int(v[12]) < 2: continue  #只保留read数目大于2的环状RNA
		circ_id = v[0] + '_' + v[1] + '_' + v[2]
		if circ_id in circ_id_lst: continue
		circ_id_lst.append(circ_id)
		f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (circ_id, v[14], v[11], v[12], v[15], v[13]))
		#circ_dict.setdefault(circ_id, {})
		#circ_dict[circ_id].setdefault(sample_id,int(v[12]))
print ('You have %d circRNAs detected' % (len(circ_id_lst)))



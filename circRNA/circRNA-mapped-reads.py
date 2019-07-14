#!/usr/bin/env python
#coding: utf-8
#从STAR mapping的LOG文件中获取uniquely mapped reads + multi mapped reads + Chimeric mapped reads

import os
import sys
#import glob
import pandas as pd
import csv


#读取shell中制作的sample name文件，按标号数字排序
sample_file = sys.argv[1]
matrix_file = sys.argv[2]


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

mapping_dict = {}
for s in files:
	sample_id = s.split('/')[1].split('@')[1]
	sample_lst.append(sample_id)
	filereader = open(s)
	line_count = 0
	total = 0
	for row in filereader:
		line_count += 1
		row = row.strip()
		item = row.split('|')
		if line_count == 9 :
			uniq_read = int(item[1].strip())
			#print(uniq_read)
			continue
		elif line_count == 24:
			multi_read = int(item[1].strip())
			#print(multi_read)
			continue
		elif line_count == 33:
			Chimeric_read = int(item[1].strip())
			#print(Chimeric_read)
			continue
		else :
			continue
	total = uniq_read + multi_read + Chimeric_read
	#print (total)
	mapping_dict.setdefault(str(sample_id), total)
#print(mapping_dict)		
df = pd.DataFrame(data = mapping_dict, index=['mapped_reads'])
#print(df)
df.to_csv(matrix_file, index=False)
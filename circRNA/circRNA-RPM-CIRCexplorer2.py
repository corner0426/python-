#!/usr/bin/env python
#coding: utf-8
#将read count矩阵数值分别除以total_mapped_reads, 得到RPM

import pandas as pd
import csv
import sys

count_file = sys.argv[1]
mapped_file = sys.argv[2]
rpm_file = sys.argv[3]

count_frame = pd.read_csv(count_file)
mapped_frame = pd.read_csv(mapped_file)

for i in count_frame.keys():
	if i in mapped_frame.keys():
		count_frame[i]=(count_frame[i]/int(mapped_frame[i]))*100000
count_frame.to_csv(rpm_file, index=False)
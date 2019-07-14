#!/bin/python
#coding:utf-8


## import modules
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd

## get mdd snp list
f1 = open('LD_clump/MDD2018_ex23.clumped', 'r')
snp_list = []
for line in f1:
    l = line.strip('\s').split()
    #print (len(l))
    if len(l) > 0 and l[0] != 'CHR':
        #snp_pos = '%s_%s' % (l[0], l[3])
        snp_id = l[2]
        #print (snp_id)
        snp_list.append(snp_id)
#print (len(snp_list)) #377

## extract mdd z statistics
mdd_A1_dict = {}
mdd_z_dict = {}
mdd_se_dict = {}
mdd_snp_pos = []
i = 0
with open('MDD2018_ex23andMe_no_sign_Zscore', 'r') as filereader:
    header = filereader.readline()
    for line in filereader:
        l = line.strip('\n').split()
        snp_pos_mdd = '%s_%s' % (l[1], l[2])
        snp_id_mdd = l[0]
        effect_allele = l[3]
        z_score = l[5]
        se = l[6]
        if snp_id_mdd in snp_list:
            i+=1
            mdd_snp_pos.append(snp_pos_mdd)
            mdd_A1_dict.setdefault(snp_pos_mdd, effect_allele)
            mdd_z_dict.setdefault(snp_pos_mdd, z_score)
            mdd_se_dict.setdefault(snp_pos_mdd, se)
print (i) # 377
mdd_dataframe = pd.DataFrame({'A1_mdd': mdd_A1_dict, 'Z_mdd':mdd_z_dict, 'se_mdd':mdd_se_dict})

## extract cur smk z statistic and make overlap file
cur_smk_A1_dict = {}
cur_smk_z_dict = {}
cur_smk_se_dict = {}
i = 0
with open('tobacco_smk_Zscore', 'r') as filereader:
    header = filereader.readline()
    for line in filereader:
        l = line.strip('\n').split()
        snp_pos = '%s_%s' % (l[1], l[2])
        snp_id = l[0]
        effect_allele = l[3]
        z_score = l[5]
        se = l[6]
        if snp_pos in mdd_snp_pos:
            i+=1
            cur_smk_A1_dict.setdefault(snp_pos, effect_allele)
            cur_smk_z_dict.setdefault(snp_pos, z_score)
            cur_smk_se_dict.setdefault(snp_pos, se)
print (i) #378
cur_smk_dataframe = pd.DataFrame({'cur_smk_A1': cur_smk_A1_dict, 'cur_smk_Z':cur_smk_z_dict, 'cur_smk_se':cur_smk_se_dict})

n = 0
l = 0
out_file = open('mdd_cur_smk_overlap.txt', 'w+')
out_file.write('CHR\tPOS\tZ_mdd\tse_mdd\tA1_mdd\tZ_cur_smk\tse_cur_smk\tA1_cur_smk\n')
for i in list(set(mdd_dataframe.index).intersection(cur_smk_dataframe.index)):
    CHR = i.split('_')[0]
    POS = i.split('_')[1]
    if mdd_dataframe['A1_mdd'][i] == cur_smk_dataframe['cur_smk_A1'][i] or mdd_dataframe['A1_mdd'][i] == Seq(cur_smk_dataframe['cur_smk_A1'][i], IUPAC.unambiguous_dna).complement():
        #(mdd_dataframe.join(cur_smk_dataframe)).to_csv(out_file)
        #print ('%s\t%s\t%s' % (mdd_dataframe['A1_mdd'][i], cur_smk_dataframe['cur_smk_A1'][i], i))
        l += 1
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], cur_smk_dataframe['cur_smk_Z'][i], cur_smk_dataframe['cur_smk_se'][i], cur_smk_dataframe['cur_smk_A1'][i]))
    else :
        n += 1
        signed_z = float(cur_smk_dataframe['cur_smk_Z'][i]) * (-1)
        out_file.write('%s\t%s\t%s\t%s\t%s\t%f\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], signed_z, cur_smk_dataframe['cur_smk_se'][i], cur_smk_dataframe['cur_smk_A1'][i]))
print (n+l) # 377

## extract past smk z statistic and make overlap file
past_smk_A1_dict = {}
past_smk_z_dict = {}
past_smk_se_dict = {}
i = 0
with open('past_tobacco_smk_Zscore', 'r') as filereader:
    header = filereader.readline()
    for line in filereader:
        l = line.strip('\n').split()
        snp_pos = '%s_%s' % (l[1], l[2])
        snp_id = l[0]
        effect_allele = l[3]
        z_score = l[5]
        se = l[6]
        if snp_pos in mdd_snp_pos:
            i+=1
            past_smk_A1_dict.setdefault(snp_pos, effect_allele)
            past_smk_z_dict.setdefault(snp_pos, z_score)
            past_smk_se_dict.setdefault(snp_pos, se)
print (i)
past_smk_dataframe = pd.DataFrame({'past_smk_A1': past_smk_A1_dict, 'past_smk_Z':past_smk_z_dict, 'past_smk_se':past_smk

n = 0
l = 0
out_file = open('mdd_past_smk_overlap.txt', 'w+')
out_file.write('CHR\tPOS\tZ_mdd\tse_mdd\tA1_mdd\tZ_past_smk\tse_past_smk\tA1_past_smk\n')
for i in list(set(mdd_dataframe.index).intersection(past_smk_dataframe.index)):
    CHR = i.split('_')[0]
    POS = i.split('_')[1]
    if mdd_dataframe['A1_mdd'][i] == past_smk_dataframe['past_smk_A1'][i] or mdd_dataframe['A1_mdd'][i] == Seq(past_smk_dataframe['past_smk_A1'][i], IUPAC.unambiguous_dna).complement():
        #(mdd_dataframe.join(cur_smk_dataframe)).to_csv(out_file)
        #print ('%s\t%s\t%s' % (mdd_dataframe['A1_mdd'][i], cur_smk_dataframe['cur_smk_A1'][i], i))
        l += 1
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], past_smk_dataframe['past_smk_Z'][i], past_smk_dataframe['past_smk_se'][i], past_smk_dataframe['past_smk_A1'][i]))
    else :
        n += 1
        signed_z = float(past_smk_dataframe['past_smk_Z'][i]) * (-1)
        out_file.write('%s\t%s\t%s\t%s\t%s\t%f\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], signed_z, past_smk_dataframe['past_smk_se'][i], past_smk_dataframe['past_smk_A1'][i]))
print (n)

## extract cur cpd z statistic and make overlap file
cur_cpd_A1_dict = {}
cur_cpd_z_dict = {}
cur_cpd_se_dict = {}
i = 0
with open('current_CPD_Zscore', 'r') as filereader:
    header = filereader.readline()
    for line in filereader:
        l = line.strip('\n').split()
        snp_pos = '%s_%s' % (l[1], l[2])
        snp_id = l[0]
        effect_allele = l[3]
        z_score = l[5]
        se = l[6]
        if snp_pos in mdd_snp_pos:
            i+=1
            cur_cpd_A1_dict.setdefault(snp_pos, effect_allele)
            cur_cpd_z_dict.setdefault(snp_pos, z_score)
            cur_cpd_se_dict.setdefault(snp_pos, se)
print (i)
cur_cpd_dataframe = pd.DataFrame({'cur_cpd_A1': cur_cpd_A1_dict, 'cur_cpd_Z':cur_cpd_z_dict, 'cur_cpd_se':cur_cpd_se_dict})

n = 0
l = 0
out_file = open('mdd_cur_cpd_overlap.txt', 'w+')
out_file.write('CHR\tPOS\tZ_mdd\tse_mdd\tA1_mdd\tZ_cur_cpd\tse_cur_cpd\tA1_cur_cpd\n')
for i in list(set(mdd_dataframe.index).intersection(cur_cpd_dataframe.index)):
    CHR = i.split('_')[0]
    POS = i.split('_')[1]
    if mdd_dataframe['A1_mdd'][i] == cur_cpd_dataframe['cur_cpd_A1'][i] or mdd_dataframe['A1_mdd'][i] == Seq(cur_cpd_dataframe['cur_cpd_A1'][i], IUPAC.unambiguous_dna).complement():
        #(mdd_dataframe.join(cur_smk_dataframe)).to_csv(out_file)
        #print ('%s\t%s\t%s' % (mdd_dataframe['A1_mdd'][i], cur_smk_dataframe['cur_smk_A1'][i], i))
        l += 1
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], cur_cpd_dataframe['cur_cpd_Z'][i], cur_cpd_dataframe['cur_cpd_se'][i], cur_cpd_dataframe['cur_cpd_A1'][i]))
    else :
        n += 1
        signed_z = float(cur_cpd_dataframe['cur_cpd_Z'][i]) * (-1)
        out_file.write('%s\t%s\t%s\t%s\t%s\t%f\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], signed_z, cur_cpd_dataframe['cur_cpd_se'][i], cur_cpd_dataframe['cur_cpd_A1'][i]))
print (n)

## extract diff not smk z statistics and make overlap file
diff_smk_A1_dict = {}
diff_smk_z_dict = {}
diff_smk_se_dict = {}
i = 0
with open('difficulty_not_smk_1day_Zscore', 'r') as filereader:
    header = filereader.readline()
    for line in filereader:
        l = line.strip('\n').split()
        snp_pos = '%s_%s' % (l[1], l[2])
        snp_id = l[0]
        effect_allele = l[3]
        z_score = l[5]
        se = l[6]
        if snp_pos in mdd_snp_pos:
            i+=1
            diff_smk_A1_dict.setdefault(snp_pos, effect_allele)
            diff_smk_z_dict.setdefault(snp_pos, z_score)
            diff_smk_se_dict.setdefault(snp_pos, se)
print (i)
diff_smk_dataframe = pd.DataFrame({'diff_smk_A1': diff_smk_A1_dict, 'diff_smk_Z':diff_smk_z_dict, 'diff_smk_se':diff_smk_se_dict})

n = 0
l = 0
out_file = open('mdd_diff_smk_overlap.txt', 'w+')
out_file.write('CHR\tPOS\tZ_mdd\tse_mdd\tA1_mdd\tZ_diff_smk\tse_diff_smk\tA1_diff_smk\n')
for i in list(set(mdd_dataframe.index).intersection(diff_smk_dataframe.index)):
    CHR = i.split('_')[0]
    POS = i.split('_')[1]
    if mdd_dataframe['A1_mdd'][i] == diff_smk_dataframe['diff_smk_A1'][i] or mdd_dataframe['A1_mdd'][i] == Seq(diff_smk_dataframe['diff_smk_A1'][i], IUPAC.unambiguous_dna).complement():
        #(mdd_dataframe.join(cur_smk_dataframe)).to_csv(out_file)
        #print ('%s\t%s\t%s' % (mdd_dataframe['A1_mdd'][i], cur_smk_dataframe['cur_smk_A1'][i], i))
        l += 1
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], diff_smk_dataframe['diff_smk_Z'][i], diff_smk_dataframe['diff_smk_se'][i], diff_smk_dataframe['diff_smk_A1'][i]))
    else :
        n += 1
        signed_z = float(diff_smk_dataframe['diff_smk_Z'][i]) * (-1)
        out_file.write('%s\t%s\t%s\t%s\t%s\t%f\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], signed_z, diff_smk_dataframe['diff_smk_se'][i], diff_smk_dataframe['diff_smk_A1'][i]))
print (n)

## extract time 1st cig z statistics and make overlap file
time_1st_A1_dict = {}
time_1st_z_dict = {}
time_1st_se_dict = {}
i = 0
with open('time_first_cig_Zscore', 'r') as filereader:
    header = filereader.readline()
    for line in filereader:
        l = line.strip('\n').split()
        snp_pos = '%s_%s' % (l[1], l[2])
        snp_id = l[0]
        effect_allele = l[3]
        z_score = l[5]
        se = l[6]
        if snp_pos in mdd_snp_pos:
            i+=1
            time_1st_A1_dict.setdefault(snp_pos, effect_allele)
            time_1st_z_dict.setdefault(snp_pos, z_score)
            time_1st_se_dict.setdefault(snp_pos, se)
print (i)
time_1st_dataframe = pd.DataFrame({'time_1st_A1': time_1st_A1_dict, 'time_1st_Z':time_1st_z_dict, 'time_1st_se':time_1st_se_dict})

n = 0
l = 0
out_file = open('mdd_time_1st_overlap.txt', 'w+')
out_file.write('CHR\tPOS\tZ_mdd\tse_mdd\tA1_mdd\tZ_time_1st\tse_time_1st\tA1_time_1st\n')
for i in list(set(mdd_dataframe.index).intersection(time_1st_dataframe.index)):
    CHR = i.split('_')[0]
    POS = i.split('_')[1]
    if mdd_dataframe['A1_mdd'][i] == time_1st_dataframe['time_1st_A1'][i] or mdd_dataframe['A1_mdd'][i] == Seq(time_1st_dataframe['time_1st_A1'][i], IUPAC.unambiguous_dna).complement():
        #(mdd_dataframe.join(cur_smk_dataframe)).to_csv(out_file)
        #print ('%s\t%s\t%s' % (mdd_dataframe['A1_mdd'][i], cur_smk_dataframe['cur_smk_A1'][i], i))
        l += 1
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], time_1st_dataframe['time_1st_Z'][i], time_1st_dataframe['time_1st_se'][i], time_1st_dataframe['time_1st_A1'][i]))
    else :
        n += 1
        signed_z = float(time_1st_dataframe['time_1st_Z'][i]) * (-1)
        out_file.write('%s\t%s\t%s\t%s\t%s\t%f\t%s\t%s\n' % (CHR, POS, mdd_dataframe['Z_mdd'][i], mdd_dataframe['se_mdd'][i], mdd_dataframe['A1_mdd'][i], signed_z, time_1st_dataframe['time_1st_se'][i], time_1st_dataframe['time_1st_A1'][i]))
print (n)



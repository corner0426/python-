#/usr/bin/python
#该脚本用于从Gendode版GTF文件中提取gene bed信息，包含基因名
import sys
import re

#path=/data/yaoyh/public_data/GTF/
gtf_file = sys.argv[1]
bed_file = sys.argv[2]

f = open(bed_file, 'w+')

with open(gtf_file, 'r') as filereader:
    for line in filereader:
        if line.startswith('#'): continue
        l = line.split('\t')
        #chr_id = l[0].strip('chr')
        start = l[3]
        end = l[4]
        if not l[0].startswith('chr'): continue 
        if l[2] == 'gene' or l[2] == 'transcript' or l[2] == 'exon':
            info = l[8].split(';')
            gene_name = 'None'
            #gene_id = 'gene'
            gene_strand = l[6]
            for i in info:
                i = i.strip(' ')
                if re.search('gene_name', i):
                    gene_name = i.split(' ')[1].strip('""')
                #elif re.search('gene_type', i):
                #    gene_status = i.split(' ')[1].strip('""')
                #elif re.search('gene_id', i):
                #    gene_id = i.split(' ')[1].strip('""').split('.')[0]
                #else:
                #    print (info)
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (l[0], start, end, l[2], gene_name, gene_strand))

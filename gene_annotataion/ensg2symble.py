#extract ensg ID and corresponding gene symble ID from gtf file 
f = open('ensg2symble.txt', 'w+')
f.write('ENSG\tgene_symble\n')

with open('Homo_sapiens.GRCh38.89_chr_version.gtf', 'r') as filereader:
    for line in filereader:
        if line.startswith('#'): continue
        ensg = ''
        symble = ''
        l = line.split('\t')
        if l[2] == 'gene':
            #print(l[8])
            ensg = l[8].split(';')[0].split(' ')[1].strip('"')
            symble = l[8].split(';')[2].split(' ')[2].strip('"')
            f.write('%s\t%s\n' % (ensg, symble))
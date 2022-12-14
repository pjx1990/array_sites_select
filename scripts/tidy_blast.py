import sys

infile = (open(sys.argv[1],'r'))
outfile = (open(sys.argv[2],'w'))

for line in infile:
    info = line.strip().split('\t')
    chrom = info[0].split(':')[0]
    qry_pos = int(info[0].split(':')[1].split('-')[0])+150
    qry_snp = chrom + '_' + str(qry_pos)
    sub_len = abs(int(info[11])-int(info[10]))
    if info[10] > info[11]:
        sub_snp = info[1]+':'+str(info[11])+'-'+str(info[10])
    else:
        sub_snp = info[1]+':'+str(info[10])+'-'+str(info[11])
    outfile.write(qry_snp+"\t"+info[0]+"\t"+sub_snp+"\t"+str(sub_len)+"\n")

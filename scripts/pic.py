import sys

weight = 40 #PIC权重

infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')

for line in infile:
    info = line.strip().split('\t')
    if info[0] == "CHROM":
        continue
    ref1 = float(info[4].split(':')[1])
    ref2 = float(info[5].split(':')[1])
    pic = 1-(ref1**2+ref2**2)
    pic_weight = int(pic * weight)
    snp = info[0]+'_'+info[1]
    outfile.write(info[0]+"\t"+info[1]+"\t"+info[2]+"\t"+info[3]+"\t"+info[4]+"\t"+info[5]+"\t"+snp+"\t"+str(pic)+"\t"+str(pic_weight)+"\n")
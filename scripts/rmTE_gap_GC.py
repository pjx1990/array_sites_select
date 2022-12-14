import sys
import re

Gap = open("1.Gap.snp",'w')
TE = open("2.TE.snp",'w')
GC = open("3.GC.snp","w")
Fa = open('snp.filter1-3.fa','w')

te_all = True  #只要有小写就过滤False，还是全部是小写时才过滤True
gc_min = 0.3
gc_max = 0.7

for line in open("snp.fa",'r'):
    info = line.strip().split('\t')
    gc = info[1].count('G')+info[1].count('C')+info[1].count('g')+info[1].count('c')
    gc_ratio = gc/len(info[1])
    chr_pos = info[0].split(":")
    pos_tmp = chr_pos[1].split("-")
    if int(pos_tmp[1]) > int(pos_tmp[0]):
        pos = int(pos_tmp[0])+150
    else:
        pos = int(pos_tmp[1])+150

    if ('N' in info[1] or 'n' in info[1]):
        continue
    else:
        Gap.write(chr_pos[0]+"\t"+str(pos)+"\n")

        if te_all is True:
            if re.search('^[a-z]+$',info[1]): #全为小写
                continue
            else:
                TE.write(chr_pos[0]+"\t"+str(pos)+"\n")
                if gc_ratio>gc_min and gc_ratio<gc_max:
                    GC.write(chr_pos[0]+"\t"+str(pos)+"\n")
                    Fa.write(">"+info[0]+"\n"+info[1]+"\n")
        else:
            if 'a' in info[1] or 'g' in info[1] or 'c' in info[1] or 't' in info[1]: #只要有小写
                continue
            else:
                TE.write(chr_pos[0]+"\t"+str(pos)+"\n")
                if gc_ratio>gc_min and gc_ratio<gc_max:
                    GC.write(chr_pos[0]+"\t"+str(pos)+"\n")
                    Fa.write(">"+info[0]+"\n"+info[1]+"\n")

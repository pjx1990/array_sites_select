import sys
import random
random.seed(0)

infile = open(sys.argv[1],'r')
print("##fileformat=VCFv4.2")
print("#CHROM"+"\t"+"POS"+"\t"+"ID"+"\t"+"REF"+"\t"+"ALT"+"\t"+"QUAL"+"\t"+"FILTER"+"\t"+"INFO"+"\t"+"FORMAT"+"\t"+"test")
for line in infile:
    info = line.strip().split("\t")
    ref = random.choice('AGCT')
    alt = random.choice('AGCT')
    print(info[0]+"\t"+info[1]+"\t"+info[2]+"\t"+ref+"\t"+alt+"\t"+str(35)+"\t"+"PASS"+"\t"+"."+"\t"+"GT"+"\t"+"1/1")

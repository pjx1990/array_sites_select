'''
GWAS候选基因。若一个基因范围内的所有变异位点数（基于背景筛选后的299894个位点注释）少于或等于5个，则全选；若多于5个，则从中选择得分最高的5个。
'''
snp_num_limit = 5

import sys
import random
import math
import heapq

random.seed(1)

snp_weight = open(sys.argv[1],'r')
ld_blocks = open(sys.argv[2],'r')
outfile = open("gwasGene.choose.snp.out",'w')

## SNP位点权重
weight = {}
for line in snp_weight:
    info = line.strip().split('\t')
    weight[info[0]] = info[1]
    # print(info[0]+":"+weight[info[0]])

## 根据字典值取键
def get_key(dct, sub_key, value):
    #字典取子集 https://www.codenong.com/3420122/
    dictfilt = lambda x, y: dict([ (i,x[i]) for i in x if i in set(y) ])
    sub_dict = dictfilt(dct, sub_key)
    k = [k for k, v in sub_dict.items() if v == value]
    # k = [k for k, v in sub_dict.items() for i in value if v == i]
    return k


## 选取得分最高的SNP
n_all = 0 #全部genes数目
n_long = 0 #超过设置阈值的genes数目
for line in ld_blocks:
    n_all += 1
    info = line.strip().split('\t')
    snps = info[1].split(' ')
    snp_len = len(snps)

    if snp_len > snp_num_limit :  
        n_long += 1
        tmp = [int(weight[key]) for key in snps] #每个gene中SNP对应的权重值
        # print(tmp)
        largest5 = heapq.nlargest(snp_num_limit,tmp) #取列表中最大的五个值，降序
        largest5.sort(reverse=True)
        # print(largest5)
        choose_id = []
        for i in largest5:
            choose_id += get_key(weight,snps,str(i))
        choose_id_rmdup = list(set(choose_id)) #去重
        ##多个SNP的分数一样时可能大于snp_num_limit 5个，则取前5个（按分数降序）。
        if len(choose_id_rmdup) > snp_num_limit:
            top5 = choose_id_rmdup[:snp_num_limit]
        else:
            top5 = choose_id_rmdup
        # print(",".join(top5))
        outfile.write(info[0]+"\t"+",".join(top5)+"\n")   
    elif snp_len <= snp_num_limit and snp_len > 0 :
        outfile.write(info[0]+"\t"+",".join(snps)+"\n")
    elif snp_len == 0:
        continue
    # break


print(f'共{n_all}个genes区间存在SNP，其中{n_long}个genes的SNP个数超过{snp_num_limit}个，从中选择得分最高的{snp_num_limit}个SNP')


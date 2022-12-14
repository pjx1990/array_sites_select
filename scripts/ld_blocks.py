'''
LD block(200kb) 计算存在LD的SNP位点，未划入任何LD block的SNP位点用于后续补缺筛选（基于均匀分布）。
如果LD block区间小于阈值（如10kb），将在其内部选择得分最高的SNP作为标签SNP，若得分最高的SNP有多个时，随机选择一个；
如果LD block大于阈值（如10kb），将在其内部每10kb（滑窗）选择一个得分最高的SNP作为这个区间的标签SNP
'''
block_len_limit = 10000

import sys
import random
import math

random.seed(123)

snp_weight = open(sys.argv[1],'r')
ld_blocks = open(sys.argv[2],'r')
ld_short = open("ld_blocks.within10kb.choose.snp.out",'w')
ld_long = open("ld_blocks.greater10kb.choose.snp.out",'w')

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
    return k

## 从LD blocks中根据距离选取得分最高的SNP
snp_num = 0  #存在LD blocks的snp总数
n_all = 0 #全部LD blocks数目
n_long = 0 #超过LD设置阈值的LD blocks数目
for line in ld_blocks:
    n_all += 1
    info = line.strip().split('\t')
    snp_num += int(info[4])
    block_len = abs(int(info[2])-int(info[1]))

    num = [int(x.split('_')[1]) for x in info[5].split('|')]
    chrx = [x.split('_')[0] for x in info[5].split('|')]
    chr_name = str(chrx[0])  #ld block只能在同一条染色体或scaffold上

    if block_len > block_len_limit: #LD block长度大于阈值
        n_long += 1
        bin_n = math.ceil(abs(num[-1]-num[0])/block_len_limit)
        # one_blocks = []
        for i in range(0,bin_n+1):
            snp_pos = [j for j in num if j >= num[0]+i*block_len_limit and j < (num[0]+(i+1)*block_len_limit)]
            # one_blocks.append(snp_pos) #空的没用，不保存了
            snp_name = [chr_name+'_'+str(i) for i in snp_pos] #加上chr/scaffold名称
            if len(snp_name) >= 1:
                tmp = [int(weight[key]) for key in snp_name] #每个ld block中SNP对应的权重值
                choose_id = get_key(weight,snp_name,str(max(tmp)))
                if len(choose_id) > 1:  #多个SNP权重最大时，从中随机取一个
                    id = random.sample(choose_id,1)
                else:
                    id = choose_id
                ld_long.write("".join(id)+'\n')

    else:  #LD block长度小于阈值
#        pass
        snp_id = info[5].split('|')
        tmp = [int(weight[key]) for key in snp_id] #每个ld block中SNP对应的权重值
        choose_id = get_key(weight,snp_id,str(max(tmp)))
        if len(choose_id) > 1:  #多个SNP权重最大时，从中随机取一个
            id = random.sample(choose_id,1)
        else:
            id = choose_id
        ld_short.write("".join(id)+'\n')
    
print(f'共{snp_num}个与其相邻SNP存在LD的SNP位点')
print(f'共{n_all}个LD blocks，其中{n_long}个LD blocks长度超过阈值，每10kb中选择一个SNP')


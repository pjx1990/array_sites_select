## softmask zs11-pan genome
bedtools maskfasta -soft -fi /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.fa -bed /mnt/project/rapeseed/data/genome/pan-zs11/zs11.v2.TE.gff -fo zs11-softmask_pan-genome.fa
makeblastdb -in zs11-softmask_pan-genome.fa -dbtype nucl -parse_seqids

## 删除重复TE、Gap和GC含量过高过低的序列
awk '{print $1"\t"$2-150"\t"$2+150}' ../02.heterozygosity/heter0.1.snp |awk '$2>0' >snp.bed
bedtools getfasta -fi /mnt/project/rapeseed/2774_filter/softmask_zs11-pan_genome/zs11-softmask_pan-genome.fa -bed snp.bed -fo snp.fa -tab
python rmTE_gap_GC.py

## 检查A、C基因组分布不均的原因
## A、C基因组分布不均由TE过滤这一步造成
#TE过滤后的SNP
sort 1.Gap.snp 2.TE.snp |uniq -u >TE.lst
#A基因组比C基因组少过滤一倍以上
grep -c 'scaffoldA' TE.lst
grep -c 'scaffoldC' TE.lst
'''
第八步筛选的SNP做均匀分布和第四步blast筛选的SNP均匀分布相比，多一些gap。分别统计这两步的gap，将第四步的SNP填补第八步。
'''

## blast比对到damor
awk -F'_' '{print $1"\t"$2-150"\t"$2+150"\t"$0}' ../all.id  |awk '$2>0' >snp.bed
bedtools getfasta -fi /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.fa -bed snp.bed -fo snp.fa
blastn -query snp.fa -db /mnt/project/rapeseed/data/genome/Darmor-bzh_v4.1/Brassica_napus_v4.1.chromosomes.fasta -out blast.out -outfmt  '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp' -num_threads 30 -max_target_seqs 10
#identity>95%, qcovs=100%, top1
awk -v OFS="\t" '$3>95 && $15=100' blast.out |awk '!a [$1]++' | awk '{print $2"_"$11"\t"$2"\t"$11}' >blast.damor4.top1.txt
cut -f 6 bed.snp.mid.txt >bed.snp.mid.lst

## snp 在damor上窗口的SNP统计
faidx /mnt/project/rapeseed/data/genome/Darmor-bzh_v4.1/Brassica_napus_v4.1.chromosomes.fasta -i chromsizes > sizes.genome
bedtools makewindows -g sizes.genome -w 100000 >damor4.100k.bed
awk '{print $2"\t"$3"\t"$1}' blast.damor4.top1.txt >chr_pos.lst
python gene_false_vcf.py chr_pos.lst >pseudo.vcf
bedtools coverage -a damor4.100k.bed -b pseudo.vcf -counts >bin.snp.count
#没有SNP的bin
awk '{if($4==0)print $0}' bin.snp.count |cut -f 1-3 >nosnp.bin
#没有SNP的bin在原始blast后SNP的bin中的SNP数
grep -f nosnp.bin /mnt/project/rapeseed/1560_filter/04.blast_rmdup/blast2damor4.1/bin.snp.count >nosnpMapBlastsnp.count
bedtools intersect -a nosnpMapBlastsnp.count2 -b /mnt/project/rapeseed/1560_filter/04.blast_rmdup/blast2damor4.1/blast.damor4.top1.bed -wa -wb >bed.snp
#用原始blast的SNP填充gap，取中间位点。且匹配回zs11-pan
Rscript get_mid.R
cat /mnt/project/rapeseed/1560_filter/04.blast_rmdup/blast2damor4.1/blast.out |awk -v OFS="\t" '$3>95 && $15=100' blast.out |awk '!a[$1]++' | awk '{print $1"\t"$2"_"$11}' >allblast.map
Rscript map.R

## 从比对回zs11-pan的region中随机取一个SNP
cat mapZS11.region |while read snp region;do
echo $region
bcftools filter /mnt/project/rapeseed/1560_filter/04.blast_rmdup/snp_biallelic.recode.vcf.gz --regions $region |grep -v '^#' |cut -f3 |shuf -n1 >>mapSNPchoose.id ##随机选取一个
done

## all SNP
cat ../all.id mapSNPchoose.id |sort -u >all.supl.gap.id


## blast damor再检验下均匀分布
awk -F'_' '{print $1"\t"$2-150"\t"$2+150"\t"$0}' ../all.supl.gap.id  |awk '$2>0' >snp.bed
bedtools getfasta -fi /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.fa -bed snp.bed -fo snp.fa
blastn -query snp.fa -db /mnt/project/rapeseed/data/genome/Darmor-bzh_v4.1/Brassica_napus_v4.1.chromosomes.fasta -out blast.out -outfmt  '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp' -num_threads 30-max_target_seqs 10
#identity>95%, qcovs=100%, top1
awk -v OFS="\t" '$3>95 && $15=100' blast.out |awk '!a [$1]++' | awk '{print $2"_"$11"\t"$2"\t"$11}' >blast.damor4.top1.txt
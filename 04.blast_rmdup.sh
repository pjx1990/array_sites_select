## blast
blastn -num_threads 40 -query snp.filter1-3.fa -db /mnt/project/rapeseed/2774_filter/softmask_zs11-pan_genome/zs11-softmask_pan-genome.fa -out blast.out -outfmt  '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp' -max_target_seqs 10

## filter: identity>90, qcovs>50
awk '{if($3>90 && $14>50)print $0}' blast.out >blast.snp.qcov.ident.filter.out
awk '{x[$1]++} END {for(i in x) print(i"\t"x[i])}' blast.snp.qcov.ident.filter.out >query_hit_num.txt
# homology hit < 5
awk '{if($2<5)print $1}' query_hit_num.txt >query_hit_num.filter
seqkit grep -f query_hit_num.filter snp.filter1-3.fa >snp.rmdup.fa
#awk '{if($2<3)print $1}' query_hit_num.txt |gzip -c >query_hit_num.filter_300K.gz

## get vcf
vcftools --gzvcf /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.id.biallelic.vcf.gz --positions blast_filter_snp_lst.txt --recode --out snp_biallelic


## zs11-pan blast2damor4.1 查看均匀分布
#grep -v '^#' ../snp_biallelic.recode.vcf |cut -f 3 >all.id
awk -F'_' '{print $1"\t"$2-150"\t"$2+150"\t"$0}' ./all.id  |awk '$2>0' >snp.bed
bedtools getfasta -fi /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.fa -bed snp.bed -fo snp.fa
blastn -query snp.fa -db /mnt/project/rapeseed/data/genome/Darmor-bzh_v4.1/Brassica_napus_v4.1.chromosomes.fasta -out blast.out -outfmt  '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp' -num_threads 30-max_target_seqs 10
#identity>95%, qcovs=100%, top1
awk -v OFS="\t" '$3>95 && $15=100' blast.out |awk '!a [$1]++' | awk '{print $2"_"$11"\t"$2"\t"$11}' >blast.damor4.top1.txt

## zs11-pan blast2damor8.1 查看均匀分布
awk -F'_' '{print $1"\t"$2-150"\t"$2+150"\t"$0}' ./all.id  |awk '$2>0' >snp.bed
bedtools getfasta -fi /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.fa -bed snp.bed -fo snp.fa
blastn -query snp.fa -db /mnt/project/rapeseed/known_array/65K_blast_damor8.1-pan/pan-Darmorv8.1.fa -out blast.out -outfmt  '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp' -num_threads 30 -max_target_seqs 10
#identity>95%, qcovs=100%, top1
awk -v OFS="\t" '$3>95 && $15=100' blast.out |awk '!a [$1]++' | awk '{print $2"_"$11"\t"$2"\t"$11}' >blast.damor8.top1.txt

## 均匀性绘图
Rscript zs11-pan_plot.R
Rscript blast2damor_plot.R


## 100kb SNP数目统计：用于后续补缺
faidx /mnt/project/rapeseed/data/genome/Darmor-bzh_v4.1/Brassica_napus_v4.1.chromosomes.fasta -i chromsizes > sizes.genome
bedtools makewindows -g sizes.genome -w 100000 >damor4.100k.bed
awk '{print $2"\t"$3"\t"$1}' blast.damor4.top1.txt >chr_pos.lst
python gene_false_vcf.py chr_pos.lst >pseudo.vcf
bedtools coverage -a damor4.100k.bed -b pseudo.vcf -counts >bin.snp.count
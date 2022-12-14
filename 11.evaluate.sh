#-------------------------------chr_distribute---------------------
## blast2damorv4.1
awk -F'_' '{print $1"\t"$2-150"\t"$2+150"\t"$0}' ../43K  |awk '$2>0' >snp.bed
bedtools getfasta -fi /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.fa -bed snp.bed -fo snp.fa
blastn -query snp.fa -db /mnt/project/rapeseed/data/genome/Darmor-bzh_v4.1/Brassica_napus_v4.1.chromosomes.fasta -out blast.out -outfmt  '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp' -num_threads 30 -max_target_seqs 10
##identity>95%, qcovs=100%, top1
#awk -v OFS="\t" '$3>95 && $15=100' blast.out |awk '!a [$1]++' | awk '{print $2"_"$11"\t"$2"\t"$11}' >blast.damor4.top1.txt
cat blast.out |awk '!a [$1]++' | awk '{print $2"_"$11"\t"$2"\t"$11}' >blast.damor4.top1.txt

## blast2damorv8.1
awk -F'_' '{print $1"\t"$2-150"\t"$2+150"\t"$0}' ../43K  |awk '$2>0' >snp.bed
bedtools getfasta -fi /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.fa -bed snp.bed -fo snp.fa
blastn -query snp.fa -db /mnt/project/rapeseed/known_array/65K_blast_damor8.1-pan/pan-Darmorv8.1.fa -out blast.out -outfmt  '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp' -num_threads 30 -max_target_seqs 10
##identity>95%, qcovs=100%, top1
#awk -v OFS="\t" '$3>95 && $15=100' blast.out |awk '!a [$1]++' | awk '{print $2"_"$11"\t"$2"\t"$11}' >blast.damor8.top1.txt
cat blast.out |awk '!a [$1]++' | awk '{print $2"_"$11"\t"$2"\t"$11}' >blast.damor8.top1.txt

## plot
Rscript zs11-pan_plot.R
Rscript blast2damor_plot.R


#-------------------------------------- LD decay -------------------------

## 43K
conda activate r4
PopLDdecay -InVCF ../variant_anno/43K.vcf -OutStat 43K
#Plot_OnePop.pl -inFile 43K.stat.gz --output 43K -bin1 10 -bin2 500
Plot_OnePop.pl -inFile 43K.stat.gz --output 43K -bin1 100 -bin2 1000

### all snp
#PopLDdecay -InVCF /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.vcf.gz -OutStat allsnp &
#Plot_OnePop.pl -inFile allsnp.stat.gz --output allsnp -bin1 100 -bin2 1000

## blast过滤后的 snp
PopLDdecay -InVCF /mnt/project/rapeseed/1560_filter/04.blast_rmdup/snp.recode.vcf -OutStat SNP
Plot_OnePop.pl -inFile SNP.stat.gz --output SNP -bin1 100 -bin2 1000


#---------------------------------MAF----------------------------------------------

vcftools --vcf /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.id.vcf --freq --out snp
vcftools --gzvcf /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_indel.vcf.gz --freq --out indel
nohup Rscript maf_stat_plot.R &


#--------------------------------PCA------------------------------------------------

conda activate r4
MingPCACluster  -InVCF ../variant_anno/43K.vcf -InSampleGroup ../ecotype_group1560.txt -OutPut 43K
Rscript 43K.plot.r
### blast后过滤的SNP
MingPCACluster  -InVCF /mnt/project/rapeseed/1560_filter/04.blast_rmdup/snp.recode.vcf -InSampleGroup ../ecotype_group1560.txt -OutPut SNP
Rscript SNP.plot.r


#--------------------------------tree----------------------------------------------

# 43K
python genetic_tree.py --geno ../variant_anno/43K.vcf --output nj_43K
# blast过滤后的SNP
python genetic_tree.py --geno /mnt/project/rapeseed/1560_filter/04.blast_rmdup/snp.recode.vcf --output nj_SNP


#--------------------------------variant annotation-----------------------------------

awk -F'_' '{print $1"\t"$2}'  ../../10.merge_allsites/43K >chr_pos.lst
bcftools view -R chr_pos.lst /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.vcf.gz >43K.snp.vcf &
awk -F'_' '{print $1"\t"$2}' ../../10.merge_allsites/43K.indel >indel.lst
bcftools view -R indel.lst /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_indel.vcf.gz >43K.indel.vcf &

grep -v '^#' 43K.indel.vcf >43K.indel.tmp  #不知为何有重复的Indel
cat 43K.snp.vcf 43K.indel.tmp >43K.vcf
java -jar snpEff/snpEff.jar -v Bnapus 43K.vcf >43K.anno.vcf

grep -v '^#' 43K.snp.vcf |cut -f 1-2,4-5 >snp.txt
grep -v '^#' 43K.indel.vcf |cut -f 1-2,4-5 |sort -u >indel.txt
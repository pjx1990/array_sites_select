##--------------------------damor sites liftover-----------------------------------------------
ref=/mnt/project/rapeseed/data/genome/Darmor-bzh_v4.1/BnIR/Brassica_napus.Darmor.v4.1.genome.fa
query=/mnt/project/rapeseed/data/genome/pan-ref/pan-reference.v0.fa
minimap2 -cx asm5 --cs $query $ref >PAF_FILE.paf
/mnt/project/Zayou_center/Denovo_project/coordinate_map/transanno/transanno-x86_64-unknown-linux-musl-v0.3.0/transanno minimap2chainPAF_FILE.paf --output zs11v0-Darv4.chain

cat ../snp.lst | cut -f 4 |sed 's/_random/-random/g' |awk -F'_' '{print $1"\t"$2"\t"$2+1"\t"$1"_"$2}' |sed 's/-random/_random/g' >damor.bed
liftOver damor.bed /mnt/project/rapeseed/data/liftOver/liftover_zs11-pan/zs11v0-Darv4.chain zs11.bed zs11_unmap.bed
cut -f 1-2 zs11.bed >zs11.pos
vcftools --gzvcf /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.vcf.gz --positions zs11.pos --recode --out zs11.snp
grep -v '^#' zs11.snp.recode.vcf |awk '{print $1"_"$2}' >snp.lst
awk '{print $1"_"$2"\t"$4}' zs11.bed >damor2zs11
grep -f snp.lst damor2zs11 >damor2zs11.snp


##--------------------------damor sites blast---------------------------------------------------

## gwas known damor sites flanking blast to zs11-pan
cat ../all.snp | sed 's/\.1//g' |sed 's/_random/-random/g' |awk -F'_'  '{print $1"\t"$2-150"\t"$2+150}' |sed 's/-random/_random/g' >damor.snp.flank150.bed
bedtools getfasta -fi /mnt/project/rapeseed/data/genome/Darmor-bzh_v4.1/BnIR/Brassica_napus.Darmor.v4.1.genome.fa -bed damor.snp.flank150.bed -fo damor.snp.flank150.fa
blastn -query damor.snp.flank150.fa -db /mnt/project/rapeseed/2774_filter/softmask_zs11-pan_genome/zs11-softmask_pan-genome.fa -out blast2pan.out -outfmt  '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp' -num_threads 30 -max_target_seqs 10 
awk '!a[$1]++{print}' blast2pan.out >blast2pan.out.top1
python tidy_blast.py blast2pan.out.top1 zs11-pan.region

## 比对zs11-pan region内的变异
#mkdir mapped
if [ -f "map.all.snp" ];then
rm map.all.snp
fi
if [ -f "damor2ZS11.snp" ];then
rm damor2ZS11.snp
fi
cat ./zs11-pan.region |while read snp darmor zs11 sub_len
do
echo $zs11
#bcftools filter /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.vcf.gz --regions $zs11 |grep -v '^#' |cut -f 1-5  >mapped/$snp
bcftools filter /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.vcf.gz --regions $zs11 |grep -v '^#' |cut -f 1-5  >>map.all.snp
##bcftools filter /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.vcf.gz --regions $zs11 |grep -v '^#' |awk '{print $1"_"$2"_"$4"_"$5}' | tr  '\n' '|' | sed 's/$/\n/' >mapped/$snp
mapsnp=`bcftools filter /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.vcf.gz --regions $zs11 |grep -v '^#' |awk '{print $1"_"$2"_"$4"_"$5}' | tr  '\n' '#' | sed 's/$/\n/'`
mapsum=`bcftools filter /mnt/project/rapeseed/data/Bna_2774/allsnp_1560_3x/1560_snp.vcf.gz --regions $zs11 |grep -v '^#' |wc -l`
echo -e "$snp\t$zs11\t$sub_len\t$mapsum\t$mapsnp" >>damor2ZS11.snp
done
sort -u map.all.snp |cut -f 1-2 >map.all.snp.uniq

## 取region内最中间的变异
Rscript gwas_choose_SNP_Indel.R
# Bnapus50K 1618 functions probes同样比对取中间变异位点
Rscript Bnapus50K_choose_SNP_Indel.R


##--------------------------genes---------------------------------------------------
## 已克隆基因
## SNP
/mnt/project/biosoft/iTools_Code/bin/iTools Gfftools AnoVar -Gff /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.gff3 -Var /mnt/project/rapeseed/1560_filter/02.heterozygosity/heter0.1.snp -OutPut snp_anno
sed 's/G/T/g' gene.lst >gene_T.lst
gunzip snp_anno.gz
grep -f gene_T.lst snp_anno |cut -f 1-2 |sort -u >uniq.snp
awk '{print $1"_"$2}' uniq.snp  >clone.gene.snp
## Indel
/mnt/project/biosoft/iTools_Code/bin/iTools Gfftools AnoVar -Gff /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.gff3 -Var /mnt/project/rapeseed/1560_filter/1560_indel_filter/chr_pos.indel -OutPut indel_anno
sed 's/G/T/g' ../gene.lst >gene_T.lst
gunzip indel_anno.gz
grep -f gene_T.lst indel_anno |cut -f 1-2 |sort -u >uniq.indel
awk '{print $1"_"$2}' uniq.indel  >clone.gene.indel

## gwas候选基因SNP
/mnt/project/biosoft/iTools_Code/bin/iTools Gfftools AnoVar -Gff /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.gff3 -Var /mnt/project/rapeseed/1560_filter/05.SNP_anno_weight/chr_pos.lst -OutPut snp_anno
gunzip snp_anno.gz
grep -f ../gene_T.lst snp_anno |cut -f 1-2 |sort -u >uniq.snp

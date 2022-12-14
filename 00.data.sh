## gatk硬过滤后的变异：2774个群体
nohup plink --vcf Bn_2774.merged.snp.filtered.vcf.gz --recode --out snp --allow-extra-chr --make-bed &

## #maf>0.05,missing rate <0.5过滤后的SNP：2774个群体
plink --bfile Bn_2774.plink.filtered  --recode vcf-iid --out Bn_2774.filtered --allow-extra-chr
bcftools view Bn_2774.filtered.vcf -Oz -o Bn_2774.filtered.vcf.gz
bcftools index Bn_2774.filtered.vcf.gz
tabix -p vcf Bn_2774.filtered.vcf.gz
plink --bfile Bn_2774.plink.filtered --recode --out Bn_2774.plink.filtered --allow-extra-chr
bcftools annotate --set-id +'%CHROM\_%POS' Bn_2774.filtered.vcf.gz >Bn_2774.filtered.id.vcf
plink --vcf Bn_2774.filtered.id.vcf --recode --out Bn_2774.plink.filtered --allow-extra-chr

## 注释
# snpEff
java -jar snpEff/snpEff.jar build -gff3 Bnapus
java -jar snpEff/snpEff.jar -v Bnapus /mnt/project/rapeseed/data/Bna_2774/filter/Bn_2774.filtered.vcf.gz >Bn_2774.anno.vcf
# iTools
grep -v '^#' ../Bn_2774.anno.vcf |cut -f 1-2 >snp.txt
/mnt/project/biosoft/iTools_Code/bin/iTools Gfftools AnoVar -Gff /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.gff3 -Var snp.txt -OutPut Bna2774 -PrintNA
/mnt/project/biosoft/iTools_Code/bin/iTools Gfftools AnoVar -Gff /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.gff3 -Var snp.txt -OutPut Bna2774_no_intergenic


## 1560个群体子集
bcftools view -S 1560_3x.lst ../Bn_2774.merged.snp.filtered.vcf.gz | bgzip -c >1560_snp.vcf.gz
bcftools view -S 1560_3x.lst ../bn_2774.merged.indel.filtered.vcf.gz | bgzip -c >1560_indel.vcf.gz
zcat 1560_snp.vcf.gz |grep -v '^#' |wc -l
#17563141

bcftools index 1560_snp.vcf.gz
tabix -p vcf 1560_snp.vcf.gz
bcftools annotate --set-id +'%CHROM\_%POS' 1560_snp.vcf.gz >1560_snp.id.vcf
plink --vcf 1560_snp.id.vcf --recode --out snp_1560 --allow-extra-chr --make-bed

bcftools view --types snps -m 2 -M 2 1560_snp.id.vcf | bgzip -c >1560_snp.id.biallelic.vcf.gz
plink --vcf 1560_snp.id.biallelic.vcf.gz --recode --out snp_1560_biallelic --allow-extra-chr --make-bed


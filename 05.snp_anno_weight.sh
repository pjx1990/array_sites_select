# 注释
ln -s /mnt/project/rapeseed/data/Bna_2774/snp_filter_maf0.05_NA0.5/annotation/snpEff .
java -jar snpEff/snpEff.jar -v Bnapus snp_biallelic.recode.vcf >snp.anno.vcf

#注释位点分级: HIGH > missense_variant > 5_prime_UTR = 3_prime_UTR > upstream_gene_variant > synonymous_variant > intron_variant > downstream_gene_variant > intergenic_region
grep "HIGH" snp.anno.vcf |awk '{print $1"_"$2"\t"20}' >HIGH.anno
cat snp.anno.vcf |grep -v "HIGH" | grep "missense_variant" |awk '{print $1"_"$2"\t"15}' >missense.anno
cat snp.anno.vcf |grep -v -E "HIGH|missense_variant" |grep -E "5_prime_UTR|3_prime_UTR" |awk '{print $1"_"$2"\t"10}' >5-3UTR.anno
cat snp.anno.vcf |grep -v -E "HIGH|missense_variant|5_prime_UTR|3_prime_UTR" |grep "upstream_gene_variant" |awk '{print $1"_"$2"\t"8}' >promoter.anno
cat snp.anno.vcf |grep -v -E "HIGH|missense_variant|5_prime_UTR|3_prime_UTR|upstream_gene_variant" |grep "synonymous_variant" |awk '{print $1"_"$2"\t"5}' >synonymous.anno
cat snp.anno.vcf |grep -v -E "HIGH|missense_variant|5_prime_UTR|3_prime_UTR|upstream_gene_variant|synonymous_variant" |grep "intron_variant" |awk '{print $1"_"$2"\t"4}' >intron.anno
cat snp.anno.vcf |grep -v -E "HIGH|missense_variant|5_prime_UTR|3_prime_UTR|upstream_gene_variant|synonymous_variant|intron_variant" | grep "downstream_gene_variant" |awk '{print $1"_"$2"\t"3}' >downstream.anno
cat snp.anno.vcf |grep -v -E "HIGH|missense_variant|5_prime_UTR|3_prime_UTR|upstream_gene_variant|synonymous_variant|intron_variant|downstream_gene_variant" |grep "intergenic_region" |awk '{print $1"_"$2"\t"2}' >intergenic.anno
#还有一些位点，基本是LOW impact，将其赋予权重8
cat snp.anno.vcf |grep -v -E "HIGH|missense_variant|5_prime_UTR|3_prime_UTR|upstream_gene_variant|synonymous_variant|intron_variant|downstream_gene_variant|intergenic_region" |grep -v '^#' |awk '{print $1"_"$2"\t"8}' >other.anno

# index
bcftools view snp_biallelic.recode.vcf  -Oz -o snp.vcf.gz
bcftools index snp.vcf.gz
cat *.anno >snp.anno.weight

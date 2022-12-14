#ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-merge_filtered_vcfs/Bn_2774.merged.snp.filtered.vcf.gz ./ -e oss-cn-beijing.aliyuncs.com -i LTAI5tSNqm6oJfQVeBdDV8YF -k rO1QCb4rQfCd7YEFuZvobNGIPBS1fT
#ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-merge_filtered_vcfs/Bn_2774.merged.snp.filtered.vcf.gz.tbi ./ -e oss-cn-beijing.aliyuncs.com -i LTAI5tSNqm6oJfQVeBdDV8YF -k rO1QCb4rQfCd7YEFuZvobNGIPBS1fT

## gatk硬过滤后的变异：2774个群体
./ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-merge_filtered_vcfs/Bn_2774.merged.snp.filtered.vcf.gz ./
./ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-merge_filtered_vcfs/Bn_2774.merged.snp.filtered.vcf.gz.tbi ./

./ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter_indel/run-1QCO3jLwtMh9uxDTKRosUln5TTe/filter_variant/5aea0811-729a-45b5-8f20-2ba1165b4360/call-merge_filtered_indel_vcfs/bn_2774.merged.indel.filtered.vcf.gz ./
./ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter_indel/run-1QCO3jLwtMh9uxDTKRosUln5TTe/filter_variant/5aea0811-729a-45b5-8f20-2ba1165b4360/call-merge_filtered_indel_vcfs/bn_2774.merged.indel.filtered.vcf.gz.tbi ./
nohup plink --vcf Bn_2774.merged.snp.filtered.vcf.gz --recode --out snp --allow-extra-chr --make-bed &

## #maf>0.05,missing rate <0.5过滤后的SNP：2774个群体
ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-merge_filtered_vcfs/Bn_2774.merged.snp.filtered.vcf.gz ./ -e oss-cn-beijing.aliyuncs.com -i LTAI5tSNqm6oJfQVeBdDV8YF -k rO1QCb4rQfCd7YEFuZvobNGIPBS1fT
ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-merge_filtered_vcfs/Bn_2774.merged.snp.filtered.vcf.gz.tbi ./ -e oss-cn-beijing.aliyuncs.com -i LTAI5tSNqm6oJfQVeBdDV8YF -k rO1QCb4rQfCd7YEFuZvobNGIPBS1fT
/mnt/project/rapeseed/data/Bna_2774/ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-plink_filter/Bn_2774.plink.filtered.bed ./
/mnt/project/rapeseed/data/Bna_2774/ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-plink_filter/Bn_2774.plink.filtered.bim ./
/mnt/project/rapeseed/data/Bna_2774/ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-plink_filter/Bn_2774.plink.filtered.fam ./
/mnt/project/rapeseed/data/Bna_2774/ossutil64 cp oss://biobin-3rdparty-saile/analysis/filter/run-1PpkbvzHlYT0zu7nJM583iR4wPR/imputation/713738f0-626b-48f2-8d58-60d4233208a6/call-plink_filter/Bn_2774.plink.filtered.log ./
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


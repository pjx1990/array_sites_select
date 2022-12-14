vcftools --gzvcf ../05.SNP_anno_weight/snp.vcf.gz --freq --out freq
python pic.py freq.frq pic.out
Rscript combine_weight.R 

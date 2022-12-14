## 有LD block的SNP根据bin window选取
cat ../../07.LD_block/all.ld_blocks.choose.snp.out |awk -F'_' '{print $1"\t"$2}' >chr_pos.id
perl ../win4snp.pl chr_pos.id chr_pos_choosed  #100kb随机取5个
awk '{print $1"_"$2}' chr_pos_choosed >ld.id

## 没有LD block的SNP根据bin window选取
grep -v '^#' ../../07.LD_block/snp.vcf |cut -f 3 >all.blast.id
cut -f 6 /mnt/project/rapeseed/1560_filter/07.LD_block/plink.blocks.det.tab |sed 's/|/\n/g' >ld_block.id
sort all.blast.id ld_block.id ld_block.id |uniq -u >nold_block.id
awk -F'_' '{print $1"\t"$2}' nold_block.id |sort -k1 -k2 -n >nold_chr_pos
vcftools --vcf /mnt/project/rapeseed/1560_filter/07.LD_block/snp.vcf  --positions nold_chr_pos --recode --out nold
grep -v '^#' nold.recode.vcf |cut -f 1-2 >nold_chr_pos2
perl win4snp_100K.pl nold.recode.vcf nold_chr_pos_choosed 100000 1  #100kb取1个
awk '{print $1"_"$2}' nold_chr_pos_choosed >nold.id

## 合并
cat ld.id nold.id >all.id

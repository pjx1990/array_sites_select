plink --bfile ../01.missing/snp.na02.maf005 --hardy --allow-extra-chr
awk '{if($7<0.1)print $2}' plink.hwe |awk -F '_' '{print $1"\t"$2}' >heter0.1.snp
/mnt/project/biosoft/iTools_Code/bin/iTools Gfftools AnoVar -Gff /mnt/project/rapeseed/data/genome/pan-zs11/pan-reference.v0.gff3 -Var heter0.1.snp -OutPut snp_anno

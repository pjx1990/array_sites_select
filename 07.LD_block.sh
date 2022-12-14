plink --vcf snp_biallelic.recode.vcf --recode --out plink --make-bed --allow-extra-chr --double-id
# plink --file plink --r2 --allow-extra-chr
# plink --bfile plink --allow-no-sex  --r2 --ld-window 999999 --ld-window-r2 0 --out test --allow-extra-chr   #https://www.jianshu.com/p/930b2743aa76
# plink --bfile plink  --show-tags test_list  --allow-extra-chr --list-all


#默认情况下仅针对200kb以内的SNP计算成对LD，可以通过--ld-window-kb选项更改此参数 #https://zhuanlan.zhihu.com/p/347399704
plink --bfile plink  --blocks no-pheno-req --allow-extra-chr
sed '1d' plink.blocks.det | sed 's/\s\+/\t/g' >plink.blocks.det.tab
#测试每个块中的每个单元型的关联性，或估算单元型的频率。需要先定相
#plink --bfile plink --hap plink.blocks --hap-freq

python ld_blocks.py ../06.PIC_weight/snp_pic.weight.txt plink.blocks.det.tab
cat ld_blocks.within10kb.choose.snp.out ld_blocks.greater10kb.choose.snp.out >all.ld_blocks.choose.snp.out

# ## A genome
# #A基因组SNP分布密度比C基因组大（重复序列过滤少）,因此LD block阈值设置高一些
# grep -E '^#|^scaffoldA' ../snp.vcf >snp_A.vcf
# plink --vcf snp_A.vcf --recode --out plink --make-bed --allow-extra-chr --double-id
# #grep '^scaffoldA' ../plink.map >plink.map
# plink --bfile plink  --blocks no-pheno-req --allow-extra-chr --ld-window-kb 400
# sed '1d' plink.blocks.det | sed 's/\s\+/\t/g' >plink.blocks.det.tab
# python ../ld_blocks.py ../../06.PIC_weight/snp_pic.weight.txt plink.blocks.det.tab
# cat ld_blocks.within10kb.choose.snp.out ld_blocks.greater10kb.choose.snp.out >all.ld_blocks.choose.snp.out
# #共172518个与其相邻SNP存在LD的SNP位点
# #共43818个LD blocks，其中1850个LD blocks长度超过阈值，每10kb中选择一个SNP
# cat all.ld_blocks.choose.snp.out ../C_genome/C.other.snp >allsnp.id

# ## C genome and other scaffold
# python ../ld_blocks.py ../../06.PIC_weight/snp_pic.weight.txt plink.blocks.det.tab
# cat ld_blocks.within10kb.choose.snp.out ld_blocks.greater10kb.choose.snp.out >all.ld_blocks.choose.snp.out
# grep -v '^scaffoldA' ../all.ld_blocks.choose.snp.out >C.other.snp
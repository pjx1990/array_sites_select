##snp：1560_background_sites.snp  50K.function.snp  clone.gene.snp  gwas.gene.snp  gwas.sites.snp
##indel：50K.function.indel  clone.gene.indel  gwas.sites.indel
cat source/1560_background_sites.snp ../09.known_sites/all.known.sites |sort -u >43K
cat source/*.snp |sort -u >43K.snp
cat source/*.indel |sort -u >43K.indel
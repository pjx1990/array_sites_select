### SNP注释和PIC权重合并
library(tidyverse)
d1 <- read.delim("pic.out",header = F)
d2 <- read.delim("../05.SNP_anno_weight/snp.anno.weight",header = F)

d <- inner_join(d2,d1[c(7,9)],by=c("V1"="V7"))
head(d)
dd <- d %>% mutate(w=V2+V9)
write.table(dd[c(1,4)],"snp_pic.weight.txt",row.names = F,col.names = F,sep = "\t",quote = F)

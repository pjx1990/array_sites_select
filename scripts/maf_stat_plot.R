#######################
## maf 统计与绘图
#######################
rm(list = ls())
library(tidyverse)


# all variant -------------------------------------------------------------


# all snp
#allsnp <- read.delim("snp.frq.test",header = F) %>% 
allsnp <- read.delim("snp.frq",header = F) %>% 
  slice(-1) %>% 
  unite(ID,c(V1,V2),sep = "_")

allsnp$V5 <- as.numeric(gsub("(.*):(.*)","\\2",allsnp$V5))
allsnp$V6 <- as.numeric(gsub("(.*):(.*)","\\2",allsnp$V6))
allsnp[allsnp=="NaN"] <- NA

allsnp <- na.omit(allsnp) %>% select(1,4,5)
allsnp$maf <- ifelse(allsnp$V5>allsnp$V6,allsnp$V6,allsnp$V5)
allsnp <- allsnp %>% select(ID,maf)

## all indel
#allindel <- read.delim("indel.frq.test",header = F) %>% 
allindel <- read.delim("indel.frq",header = F) %>% 
  slice(-1) %>% 
  unite(ID,c(V1,V2),sep = "_")

allindel$V5 <- as.numeric(gsub("(.*):(.*)","\\2",allindel$V5))
allindel$V6 <- as.numeric(gsub("(.*):(.*)","\\2",allindel$V6))
allindel[allindel=="NaN"] <- NA

allindel <- na.omit(allindel) %>% select(1,4,5)
allindel$maf <- ifelse(allindel$V5>allindel$V6,allindel$V6,allindel$V5)
allindel <- allindel %>% select(ID,maf)

## all
allsnp$Type <- "SNP"
allindel$Type <- "Indel"
allvariant <- rbind(allsnp,allindel)
write.table(allvariant,"allvariant.maf",row.names = F,col.names = T,sep = "\t",quote = F)

## 统计绘图
break_maf <- c(0,0.05,0.1,0.2,0.3,0.4,0.5)
labels = c("[0,0.05]", "(0.05,0.1]", "(0.1,0.2]", "(0.2,0.3]","(0.3,0.4]","(0.4,0.5]")
pd <- data.frame(table(cut(allvariant$maf,break_maf,labels)))
names(pd) <- c("MAF","Count")
p=ggplot(pd,aes(MAF,Count,fill=Count))+
  geom_bar(stat = "identity")+
  geom_col(position = "dodge")+
  geom_text(aes(label = Count), position = position_dodge(0.5), vjust = -0.5)+
  theme(legend.position = "none");p
ggsave(p,filename = "allvariant.png",width = 6,height = 4,dpi = 300)
ggsave(p,filename = "allvariant.pdf",width = 6,height = 4,dpi = 300)


# 43K ---------------------------------------------------------------------

## snp
snp43K <- read.delim("../../10.merge_allsites/43K.snp",header = F)
snpmap <- inner_join(allsnp,snp43K,by=c("ID"="V1"))

## indel
indel43K <- read.delim("../../10.merge_allsites/43K.indel",header = F)
indelmap <- inner_join(allindel,indel43K,by=c("ID"="V1"))

## all
variant43K <- rbind(snpmap,indelmap)
write.table(variant43K,"variant43K.maf",row.names = F,col.names = T,sep = "\t",quote = F)

## 统计绘图
break_maf <- c(0,0.05,0.1,0.2,0.3,0.4,0.5)
labels = c("[0,0.05]", "(0.05,0.1]", "(0.1,0.2]", "(0.2,0.3]","(0.3,0.4]","(0.4,0.5]")
pd <- data.frame(table(cut(variant43K$maf,break_maf,labels)))
names(pd) <- c("MAF","Count")
p=ggplot(pd,aes(MAF,Count,fill=Count))+
  geom_bar(stat = "identity")+
  geom_col(position = "dodge")+
  geom_text(aes(label = Count), position = position_dodge(0.5), vjust = -0.5)+
  theme(legend.position = "none");p
ggsave(p,filename = "variant43K.png",width = 6,height = 4,dpi = 300)
ggsave(p,filename = "variant43K.pdf",width = 6,height = 4,dpi = 300)

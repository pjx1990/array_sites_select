#########################
## 筛选damor blast 后的ZS11-pan中的位点（SNP/Indel）
## 根据zs11靠中间位置，且ref alt与damor位点相同
#########################
rm(list = ls())
library(tidyverse)

## liftover SNP and Indel
snp <- read.delim("damor2zs11.snp",header = F)
indel <- read.delim("damor2zs11.indel",header = F)
snp$Target_Type <- "SNP"
indel$Target_Type <- "INDEL"
variant <- rbind(snp,indel)
names(variant)[1:2] <- c("ZS11_ID","Damor_ID")

##固相50K所有功能探针
d1 <- read.csv("1619function_probes.csv",header = T) %>% 
  filter(Chromosome != "other") %>% 
  filter(Chromosome != "mtchr") %>% 
  select(4,5,6,8)
names(d1)

##加上变异类型和ref alt测试下
liftv <- inner_join(d1,variant,by=c("Chr_pos"="Damor_ID","Target_Type"))
write.table(liftv,"liftover_success.txt",sep = "\t",quote = F,row.names = F,col.names = T)

## 剩余的SNP取blast比对的中间位点
left_snp <- setdiff(d1$Chr_pos,variant$Damor_ID)
dd1 <- d1[match(left_snp,d1$Chr_pos),]

##比对结果
d2 <- read.delim("damorFlanking150_blast_ZS11pan.snp.indel",header = F)

d3 <- inner_join(dd1,d2,by=c("Chr_pos"="V1"))
names(d3) <- c("Damor_ID","Damor_Ref","Damor_Alt","Variant_Type","ZS11_Region","Subject_Len","Variant_num","Variant")

set.seed(0)
choose_snp <- NULL
for(i in 1:nrow(d3)){
  # i=1
  # i=4
  snps <- unlist(strsplit(d3[i,8],"#"))
  if(d3[i,7]==0){
    snp <- NA
  }else if(d3[i,7]%%2==0){ #偶数个SNP，选择中间两个（优先选择与damor ref和alt一致的位点，否则随机选一个）。##否则选第一个
    ix1 <- d3[i,7]/2
    snp1 <- snps[ix1]
    snp_info1 <- unlist(strsplit(snp1,"_"))
    ix2 <- (d3[i,7]/2)+1
    snp2 <- snps[ix2]
    snp_info2 <- unlist(strsplit(snp2,"_"))
    if(d3[i,2]==snp_info1[3] & d3[i,3]==snp_info1[4]){
      snp <- snp1
    }else if(d3[i,2]==snp_info1[3] & d3[i,3]==snp_info2[4]){  
      snp <- snp2
    }else{
      ix <- sample(c(ix1,ix2),1)
      # ix <- ix1
      snp <- snps[ix]
    }
    
  }else if(d3[i,7]%%2==1){ #奇数个SNP
    ix <- (d3[i,7]+1)/2
    snp <- snps[ix]
  }
  choose_snp <- c(choose_snp,snp)
}

d4 <- d3
d4$Choosed_Variant <- choose_snp
d5 <- d4 %>% separate(Choosed_Variant,c("ZS11_Chr","ZS11_Pos","ZS11_Ref","ZS11_Alt"),sep = "_") %>% 
  unite(Chr_Pos,c(ZS11_Chr,ZS11_Pos),sep = "_")
d5$Chr_Pos <- gsub("NA_NA",NA,d5$Chr_Pos)
d5$Choosed_Variant <- choose_snp

write.csv(d5,"damor2ZS11_choosed_variant.csv",row.names = F)


# all variant list --------------------------------------------------------

snpd1 <- liftv %>% filter(Target_Type=="SNP")
snpd2 <- d5 %>% filter(Variant_Type=="SNP")
snp_list <- unique(c(snpd1$ZS11_ID,snpd2$Chr_Pos)) 
snp_list <- snp_list[!is.na(snp_list)]
write.table(snp_list,"all_uniq_SNP_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)


indeld1 <- liftv %>% filter(Target_Type=="INDEL")
indeld2 <- d5 %>% filter(Variant_Type=="INDEL")
indel_list <- unique(c(indeld1$ZS11_ID,indeld2$Chr_Pos)) 
indel_list <- indel_list[!is.na(indel_list)]
write.table(indel_list,"all_uniq_Indel_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)

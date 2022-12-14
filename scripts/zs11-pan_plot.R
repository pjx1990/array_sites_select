library(CMplot)
library(tidyverse)

# bin -----------------------------------------------------------------


data3 <-  read.delim("43K",header = F) 
dd <- data3 %>% separate(V1,c("V2","V3"),sep = "_")
dd2 <- cbind(data3,dd) %>% filter(!grepl("scaffold[0-9]",V2)) %>% arrange(V2)
names(dd2) <- c("SNP", "Chromosome", "Position")
head(dd2)
# table(data3$Chromosome)



# 1M window ---------------------------------------------------------------


CMplot(
  dd2, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="window1M", dpi=300, file.output=TRUE, verbose=TRUE
)



dd3 <- dd2 %>% filter(!grepl("genome|damor",Chromosome))
CMplot(
  dd3, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="window1M_scaffold", dpi=300, file.output=TRUE, verbose=TRUE
)



# 100k window ---------------------------------------------------------------


CMplot(
  dd2, plot.type="d",  bin.size=1e5, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="window100k", dpi=300, file.output=TRUE, verbose=TRUE
)



dd3 <- dd2 %>% filter(!grepl("genome|damor",Chromosome))
CMplot(
  dd3, plot.type="d",  bin.size=1e5, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="_scaffold_window100k", dpi=300, file.output=TRUE, verbose=TRUE
)


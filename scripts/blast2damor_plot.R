library(CMplot)
library(tidyverse)

# 43K blast damor v4.1---------------------------------------------------------------------


data <- read.delim("blast.damor4.top1.txt",header = F) %>% arrange(V2)
names(data) <- c("SNP", "Chromosome", "Position")
unique(data$Chromosome)

CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="balst2damor4", dpi=300, file.output=TRUE, verbose=TRUE
)

d2 <- data %>% filter(!grepl("random",Chromosome))
CMplot(
  d2, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="balst2damor4_chr", dpi=300, file.output=TRUE, verbose=TRUE
)

CMplot(
  data, plot.type="d",  bin.size=1e5, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="balst2damor4_win100K", dpi=300, file.output=TRUE, verbose=TRUE
)

d2 <- data %>% filter(!grepl("random",Chromosome))
CMplot(
  d2, plot.type="d",  bin.size=1e5, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="balst2damor4_chr_win100K", dpi=300, file.output=TRUE, verbose=TRUE
)
CMplot(
  d2, plot.type="d",  bin.size=1e4, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="balst2damor4_chr_win10K", dpi=300, file.output=TRUE, verbose=TRUE
)


# 43K blast damor v8.1 pan---------------------------------------------------------------------


data <- read.delim("blast.damor8.top1.txt",header = F) %>% arrange(V2)
names(data) <- c("SNP", "Chromosome", "Position")
unique(data$Chromosome)

CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="blast2damor8Pan", dpi=300, file.output=TRUE, verbose=TRUE
)

d2 <- data %>% filter(grepl("chr",Chromosome))
CMplot(
  d2, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="blast2damor8Pan_chr", dpi=300, file.output=TRUE, verbose=TRUE
)


CMplot(
  data, plot.type="d",  bin.size=1e5, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="blast2damor8Pan_win100K", dpi=300, file.output=TRUE, verbose=TRUE
)

d2 <- data %>% filter(grepl("chr",Chromosome))
CMplot(
  d2, plot.type="d",  bin.size=1e5, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="blast2damor8Pan_chr_win100K", dpi=300, file.output=TRUE, verbose=TRUE
)

CMplot(
  d2, plot.type="d",  bin.size=1e4, col=c("darkgreen", "yellow", "red"),
  file="jpg", memo="blast2damor8Pan_chr_win10K", dpi=300, file.output=TRUE, verbose=TRUE
)

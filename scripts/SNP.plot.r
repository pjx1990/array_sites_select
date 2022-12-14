if(!require(gplots)) install.packages("gplots")
if(!require(ggplot2)) install.packages("ggplot2")
library(gplots)
library(ggplot2)
read.table("SNP.eigenvec",h=T)->PCA

grc = factor(PCA[["Group"]])
shape_level <- nlevels(grc)
if (shape_level < 15){ shapes = (0:shape_level) %% 15}else{ shapes = c(0:14,c((15:shape_level) %% 110 + 18)) }


pdf("SNP.r.PCA1_PCA2.pdf")
if (shape_level < 6 ){
ggplot(data = PCA, mapping = aes(x = PCA1, y = PCA2,colour=Cluster,shape =Group)) + geom_point(size = 3 ) +ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5)) 
} else {
ggplot(data = PCA, mapping = aes(x = PCA1, y = PCA2,colour=Cluster,shape =Group)) + geom_point(size = 3 ) +ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5))+  scale_shape_manual(values=seq(0,15)) 	}

dev.off();

pdf("SNP.r.PCA1_PCA3.pdf")
if (shape_level < 6 ){
ggplot(data = PCA, mapping = aes(x = PCA1, y = PCA3,colour=Cluster,shape =Group)) + geom_point(size = 3) +ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5)) 
} else {
ggplot(data = PCA, mapping = aes(x = PCA1, y = PCA3,colour=Cluster,shape =Group)) + geom_point(size = 3) +ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5))  +  scale_shape_manual(values=seq(0,15))
}
dev.off();


read.table("SNP.evaluation",h=T)->SSE
pdf("SNP.K_SSE.pdf")
ggplot(data = SSE, mapping = aes(x =K,y=SSE),color=season) +geom_line(size=2,linetype=1)+geom_point(size = 3) 
dev.off()
pdf("SNP.K_DBi.pdf")
ggplot(data = SSE, mapping = aes(x =K,y=DBi),color=season) +geom_line(size=2,linetype=1)+geom_point(size = 3)
dev.off()

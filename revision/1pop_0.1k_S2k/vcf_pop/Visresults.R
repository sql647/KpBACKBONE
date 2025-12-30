library(ggtree)
library(treeio)
library(cowplot)
library(ggplot2)

t<-read.newick("pop1_5000.FastTree")
ggtree(t,layout = "ape")+
  geom_treescale(offset = 0.005,fontsize=5, linesize=1, x=0, y=0,width=0.01)

gd<-read.table("../gamma_distribution_effect.txt")
colnames(gd)<-"coefficien"
hist(gd$coefficien,xlab = "Coeffcient",main = "Histogram")
gd$POSITION<-c(1:length(gd$coefficien))

plot(gd$POSITION,gd$coefficien,type = "p",
     xlab="position",ylab="coefficent",main = "rgamma(1000000,shape = 0.01,scale = 100)")

###vis_PCA
pca <- read.table("pop1maf002pca.eigenvec")
eigenval <- scan("pop1maf002pca.eigenval")
# remove nuisance column
rownames(pca)<-pca$V1
pca <- pca[,-c(1,2)]
colnames(pca) <- paste0("PC", 1:ncol(pca))
val<-read.table("pop1maf002pca.eigenval")
colnames(val)<-"val"
rownames(val)<-paste0("PC", 1:nrow(val))
pca<-as.data.frame(t(t(pca)*sqrt(val$val)))
# set names
pca$ID<-rownames(pca)
##pve
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
library(ggplot2)
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
  ylab("Percentage variance explained") + theme_light()+
  theme( 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
    axis.title.y = element_text(size = 16),   # 设置y轴标题的字体大小（如果需要的话）
    plot.title = element_text(size = 20)
  )+
  xlab("Principle Component")+
  ggtitle("Scree Plot")
###
###
ggplot(pca,aes(x=PC1,y=PC2))+
  geom_point(size=5)+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.size = unit(2, "lines"))
###
library(grid)
a<-ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
  ylab("PVE") + theme_light()+
  theme( 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 12),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 12),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 13),  # 设置x轴标题的字体大小（如果需要的话）  
    axis.title.y = element_text(size = 13),   # 设置y轴标题的字体大小（如果需要的话）
    plot.title = element_text(size = 12)
  )+
  xlab("Principal Component")+
  ggtitle("Scree Plot")

scree_grob <- ggplotGrob(a)

ggplot(pca,aes(x=PC1,y=PC2))+
  geom_point(size=5,alpha = 0.5)+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.size = unit(2, "lines"))+
  annotation_custom(grob = scree_grob, xmin = -0.15, xmax = 1.85, ymin = -1, ymax = 1)

###snp
geno<-read.table("pop1maf002.matrix.MajorAMinorT",header = TRUE, row.names = NULL)
pos <- geno[, 1, drop = FALSE]  # 使用[, 1]选择第一列，drop = FALSE保持为数据框  
geno <- geno[, -1]  # 移除第一列（位置信息），只保留基因型数据  
geno[geno=="A"]<- 0
geno[geno=="T"]<- 1
geno <- data.frame(lapply(geno, as.numeric))
load_snp <- load_strain %*% t(geno)
load_snp_squ<-load_snp^2

load_snp_df<-as.data.frame(t(load_snp))
load_snp_df$POS<-as.numeric(pos$row.names)

squdf<-as.data.frame(t(load_snp_squ))
squdf$POS<-as.numeric(pos$row.names)

plot(load_snp_df$POS,load_snp_df$PC1,type = "p",xlab = "POSITION",ylab = "PC1 SNP Loading")

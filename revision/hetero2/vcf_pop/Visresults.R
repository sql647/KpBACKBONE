library(ggtree)
library(treeio)
library(cowplot)
library(ggplot2)

t<-read.newick("pop1_10000.FastTree")
ggtree(t,layout = "ape")+
  geom_treescale(offset = 0.005,fontsize=5, linesize=1, x=0, y=0,width=0.01)


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
#
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
  geom_point(size=5,alpha=0.5)+
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
  annotation_custom(grob = scree_grob, xmin = -1, xmax = -0.2, ymin = -0.9, ymax = -0.1)

###
##strain loadings
pca<-read.table("pop1maf002pca.eigenvec")
rownames(pca)<-pca$V1
pca$V1<-NULL
pca$V2<-NULL
colnames(pca)<-paste0("PC",1:20)
val<-read.table("pop1maf002pca.eigenval")
colnames(val)<-"val"
rownames(val)<-paste0("PC", 1:nrow(val))
load_strain<-t(pca)*sqrt(val$val)
###snp
geno<-read.table("pop1maf002.matrix",header = TRUE, row.names = NULL)
pos <- geno[, 1, drop = FALSE]
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
plot(load_snp_df$POS,load_snp_df$PC1,type = "p")

ggplot(squdf,aes(x=POS,y=PC1))+
  geom_point( alpha = 0.5, size = 3)+
  labs(x = "SNP Position", y = "Squared PC1 Loading") +
  theme_classic()+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")+
  geom_vline(xintercept = 250000,color="red")+geom_vline(xintercept = 270000,color="red")+
  geom_vline(xintercept = 500000,color="red")+geom_vline(xintercept = 520000,color="red")+
  geom_vline(xintercept = 750000,color="red")+geom_vline(xintercept = 770000,color="red")
###
###
library("CMplot")

cmdata <- data.frame(
  SNP = paste0("ID", seq_len(nrow(squdf))),
  CHR = 1,                  # 所有点都属于同一条染色体 → 不会丢点
  BP  = squdf$POS,          # 真实 bp 位置，一点不丢
  P   = squdf$PC1           # loading，全部保留
)

CMplot(
  cmdata,
  plot.type = "c",
  type = "p",
  LOG10 = FALSE,
  r = 2,     # 控制中间空心大小
  col = "#00000080",
  signal.line = NULL,
  outward = TRUE,
  cex = 1,
  signal.cex = 1,
  cir.chr.h=.5,
  chr.labels = FALSE,  # 不让 CMplot 自动写染色体号（因为只有 CHR=1）
  dpi = 600,
  file = "pdf",
  file.output = TRUE
)
  

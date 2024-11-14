###vis_PCA
pca <- read.table("NR369pca.eigenvec")
eigenval <- scan("NR369pca.eigenval")
# remove nuisance column
rownames(pca)<-pca$V1
pca <- pca[,-c(1,2)]
colnames(pca) <- paste0("PC", 1:ncol(pca))
val<-read.table("NR369pca.eigenval")
colnames(val)<-"val"
rownames(val)<-paste0("PC", 1:nrow(val))
pca<-as.data.frame(t(t(pca)*sqrt(val$val)))
# set names
pca$ID<-rownames(pca)
##pve
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
library(ggplot2)
library(ggsci)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
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


library(grid)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
  ylab("PVE") + theme_light()+
  theme( 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 7),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 7),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 8),  # 设置x轴标题的字体大小（如果需要的话）  
    axis.title.y = element_text(size = 8),   # 设置y轴标题的字体大小（如果需要的话）
    plot.title = element_text(size = 10)
  )+
  xlab("Principle Component")+
  ggtitle("Scree Plot")
scree_grob <- ggplotGrob(a)
ggplot(pca_info,aes(x=PC1,y=PC2))+
  geom_point(aes(size=virulence_score, 
                 shape=factor(resistance_score),
                 color=Continent))+
  scale_color_aaas()+
  theme_light()+
  scale_shape_discrete(name = "Resistance Score")+
  scale_size_continuous(name = "Virulence Score",transform = "reverse",range = c(1, 3))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE,title.position = "top"),
         size=guide_legend(nrow=2, byrow=TRUE,title.position = "top"),
         shape=guide_legend(nrow=2, byrow=TRUE,title.position = "top"))+
  theme(legend.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "lines"),
        legend.position = "bottom")+
  annotation_custom(grob = scree_grob, xmin = -0.1, xmax = 0.2, ymin = 0.1, ymax = 0.4)

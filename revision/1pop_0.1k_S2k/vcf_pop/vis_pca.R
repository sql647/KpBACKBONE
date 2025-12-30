library(tidyverse)
pca <- read_table("pop1maf002pca.eigenvec", col_names = FALSE)
eigenval <- scan("pop1maf002pca.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
##pve
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

##make plot
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

# plot pca
b <- ggplot(pca, aes(PC1, PC2,color=PC1)) + 
  geom_point(size = 3)+
  theme_light()+
  scale_color_gradientn(colors = c("#FF0000" ,"#0000FF","#A0E8AF","#000000"))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16),
        legend.position = "none")

library(cowplot)
plot_grid(b,a, align = "hv")

###
library(ggtree)
library(treeio)
pop1tree<-read.newick("pop1_10000.FastTree")
colnames(pca)[1]<-"label"

tree.a<-full_join(pop1tree,pca,by="label")

ggtree(tree.a, layout = "ape") +  
  geom_tippoint(aes(color = PC1),size=3) +  
  scale_color_gradientn(  
    colours = c("#FF0000" ,"#0000FF","#A0E8AF","#000000"), # 高对比度颜色  
    values = scales::rescale(seq(min(pca$PC1), max(pca$PC1), length.out = 5)), # 根据PC1的值调整颜色对应值  
    name = "PC1" # 图例标题  
  )+
  geom_treescale(offset = 0.01,fontsize=5, linesize=1, x=0, y=0,width=0.01)+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(2, "lines"),
        legend.position = "bottom")

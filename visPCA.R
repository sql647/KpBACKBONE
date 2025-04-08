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

metainfo<-read.csv("allKpmeta.csv")
library(dplyr)
pca_info <-left_join(pca, metainfo)
countryinfo <- strsplit(as.character(pca_info$geo_loc_name), ":")
countryinfo1 <- sapply(countryinfo, `[`, 1)
pca_info$country<-countryinfo1
library(countrycode)
pca_info$Continent<-countrycode(pca_info$country, "country.name", "continent")
pca_info$Continent <- ifelse(is.na(pca_info$Continent), "NA", pca_info$Continent)

###
y<-pca_info[,c("ID","PC1","Continent","virulence_score","resistance_score",
               "geo_loc_name","host","isolation_source","host_disease","host_health_state")]
y[297:298,"Continent"]<-"Europe"
y$Continent <- ifelse(is.na(y$Continent), "missing", y$Continent)
y$Continent <- ifelse(is.na(y$Continent), "missing", y$Continent)
y$host[y$host == ""] <- "missing"
y$host[y$host == "not collected"] <- "missing"
y$host[y$host == "Not Applicable"] <- "missing"
y$host[y$host == "not available"] <- "missing"
y$host[y$host == "not applicable" ] <- "missing"
y$`geo_loc_name`[y$`geo_loc_name` == ""] <- "missing"
y$ID <- sub("_(\\d+)$", ".\\1", y$ID)
y<-y[order(y$PC1),]
#write.csv(y,file = "figs/NR369info.csv",quote = FALSE,row.names = FALSE)

x<-pca_info[,c("ID",paste0("PC",1:20),"Continent","virulence_score","resistance_score")]
x$Continent <- ifelse(is.na(x$Continent), "missing", x$Continent)
#write.csv(x,file = "fs/NR369info.csv",quote = FALSE,row.names = FALSE)
library("scales")
show_col(pal_aaas("default")(10))
b <- ggplot(pca_info,aes(x=PC1,y=PC2,color=Continent))+
  geom_point(size=5)+
  theme_light()+
  scale_color_aaas()+
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
b
###add info
pca_info[297:298,"Continent"]<-"Europe"
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
        legend.position = "bottom")
##host
library(stringr)

pca_info$host_clean <- tolower(trimws(pca_info$host))

pca_info$host_clean <- case_when(
  str_detect(pca_info$host_clean, "homo|human") ~ "Homo sapiens",
  str_detect(pca_info$host_clean, "swine|pig|boar") ~ "Sus scrofa",
  str_detect(pca_info$host_clean, "cow|cattle|calf|bos taurus") ~ "Bos taurus",
  str_detect(pca_info$host_clean, "dog|canis lupus familiaris") ~ "Canis lupus familiaris",
  str_detect(pca_info$host_clean, "cat|felis catus") ~ "Felis catus",
  str_detect(pca_info$host_clean, "mouse|mus musculus") ~ "Mus musculus",
  str_detect(pca_info$host_clean, "chicken|gallus gallus") ~ "Gallus gallus",
  str_detect(pca_info$host_clean, "monkey") ~ "Monkey",
  str_detect(pca_info$host_clean, "hedgehog") ~ "Hedgehog",
  str_detect(pca_info$host_clean, "sheep|ovis aries") ~ "Ovis aries",
  str_detect(pca_info$host_clean, "horse|equus caballus") ~ "Equus caballus",
  str_detect(pca_info$host_clean, "bat|pteropus") ~ "Pteropus poliocephalus",
  str_detect(pca_info$host_clean, "phocarctos hookeri") ~ "Phocarctos hookeri",
  str_detect(pca_info$host_clean, "garcinia xanthochymus|plant") ~ "Plant",
  str_detect(pca_info$host_clean, "environment|water|housefly") ~ "Environment",
  str_detect(pca_info$host_clean, "not|unknown|missing|^$") ~ "Unknown",
  TRUE ~ "Unknown"
)

pca_info$host_simple <- dplyr::recode(pca_info$host_clean,
                                      "Homo sapiens" = "Human",
                                      "Sus scrofa" = "Pig",
                                      "Bos taurus" = "Cow",
                                      "Mus musculus" = "Mouse",
                                      "Gallus gallus" = "Chicken",
                                      "Monkey" = "Monkey",
                                      "Ovis aries" = "Sheep",
                                      "Equus caballus" = "Horse",
                                      "Canis lupus familiaris" = "Dog",
                                      "Felis catus" = "Cat",
                                      "Phocarctos hookeri" = "Sea lion",
                                      "Pteropus poliocephalus" = "Bat",
                                      "Hedgehog" = "Hedgehog",
                                      "Plant" = "Plant",
                                      "Environment" = "Environment",
                                      "Unknown" = "Unknown"
)
###
# 生成 16 种彩虹色
host_colors <- rainbow(n = 16)

# 给颜色命名（确保顺序一致）
names(host_colors) <- sort(unique(pca_info$host_simple))

###
library(grid)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
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
ggplot(pca_info,aes(x=PC1,y=PC2))+
  geom_point(aes(size=virulence_score, 
                 shape=factor(resistance_score),
                 color=Continent,
                 fill=host_simple))+
  scale_color_npg()+
  theme_light()+
  scale_shape_manual(
    name = "Resistance Score",
    values = c("0" = 22, "1" = 23, "2" = 24, "3" = 25)
  )+
  scale_size_continuous(name = "Virulence Score",transform = "reverse",range = c(1, 3))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  #ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  guides(color=guide_legend(nrow=3, byrow=TRUE,title.position = "top"),
         size=guide_legend(nrow=2, byrow=TRUE,title.position = "top"),
         shape=guide_legend(nrow=1, byrow=TRUE,title.position = "top"),
         fill = guide_legend(
           title.position = "top",
           nrow = 8,
           byrow = TRUE,
           override.aes = list(shape = 21, color = "black", stroke = 1)
         ))+
  scale_fill_manual(name = "Host", values = host_colors)+
  theme(legend.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"),
        legend.position = "right")+
  annotation_custom(grob = scree_grob, xmin = -0.07, xmax = 0.23, ymin = 0.1, ymax = 0.4)+
  scale_x_continuous(breaks = seq(-0.3, 0.4, by = 0.2))+
  geom_density(aes(y = ..density.. /10), fill = "skyblue", alpha = 0.2,color = NA) +  
  scale_y_continuous(
    name = paste0("PC2 (", signif(pve$pve[2], 3), "%)"),
    sec.axis = sec_axis(~ . * 10, name = "Density of PC1")  
  )

ggplot(pca_info,aes(x=PC1,y=PC2))+
  geom_point(aes(size=virulence_score, 
                 shape=factor(resistance_score),
                 color=Continent))+
  scale_color_npg()+
  theme_light()+
  scale_shape_discrete(name = "Resistance Score")+
  scale_size_continuous(name = "Virulence Score",transform = "reverse",range = c(1, 3))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  #ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  guides(color=guide_legend(nrow=3, byrow=TRUE,title.position = "top"),
         size=guide_legend(nrow=2, byrow=TRUE,title.position = "top"),
         shape=guide_legend(nrow=1, byrow=TRUE,title.position = "top"))+
  theme(legend.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "lines"),
        legend.position = "right")+
  annotation_custom(grob = scree_grob, xmin = -0.07, xmax = 0.23, ymin = 0.1, ymax = 0.4)+
  scale_x_continuous(breaks = seq(-0.3, 0.4, by = 0.2))+
  geom_density(aes(y = ..density.. /10), fill = "skyblue", alpha = 0.2,color = NA) +  
  scale_y_continuous(
    name = paste0("PC2 (", signif(pve$pve[2], 3), "%)"),
    sec.axis = sec_axis(~ . * 10, name = "Density of PC1")  
  )

###
library(GGally)
pc.pr <- ggpairs(pca_info,columns = c(1:5),
                 mapping = ggplot2::aes(color=Continent),
                 upper = list(continuous = "points"),
                 diag = list(continuous = "blankDiag"))+
  scale_color_aaas()+
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        panel.background = element_blank(),)
pc.pr

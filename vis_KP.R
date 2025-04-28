
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
library(writexl)

#write_xlsx(y, path = "../peer review/suptable1.xlsx")

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
    axis.title.x = element_text(size = 14),  # 设置x轴标题的字体大小（如果需要的话）  
    axis.title.y = element_text(size = 14),   # 设置y轴标题的字体大小（如果需要的话）
    plot.title = element_text(size = 14)
  )+
  xlab("PC")+
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

fig1b<-ggplot(pca_info,aes(x=PC1,y=PC2))+
  geom_point(aes(size=virulence_score, 
                 shape=factor(resistance_score),
                 color=Continent))+
  scale_color_npg()+
  theme_light()+
  scale_shape_discrete(name = "Resistance Score")+
  scale_size_continuous(name = "Virulence Score",transform = "reverse",range = c(1, 3))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  #ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
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
  annotation_custom(grob = scree_grob, xmin = -0.07, xmax = 0.23, ymin = 0.1, ymax = 0.4)+
  scale_x_continuous(breaks = seq(-0.3, 0.4, by = 0.2))+
  geom_density(aes(y = ..density.. /10), fill = "skyblue", alpha = 0.2,color = NA) +  
  scale_y_continuous(
    name = paste0("PC2 (", signif(pve$pve[2], 3), "%)"),
    sec.axis = sec_axis(~ . * 10, name = "Density of PC1")  
  )
ggsave("/Users/47liu/Documents/paper/peer review/figs/figure1b.pdf",
       plot = fig1b,
       width = 8, height = 6,
       units = "in",
       dpi = 600,
       bg = "transparent",
       useDingbats = FALSE  
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

###loadings
pca<-read.table("NR369pca.eigenvec")
rownames(pca)<-pca$V1
pca$V1<-NULL
pca$V2<-NULL
colnames(pca)<-paste0("PC",1:20)
val<-read.table("NR369pca.eigenval")
colnames(val)<-"val"
rownames(val)<-paste0("PC", 1:nrow(val))
##strain loadings
load_strain<-t(pca)*sqrt(val$val)

###gene
genetab<-read.table("NRdata_gene_maf002.Rtab",check.names = FALSE,header = TRUE)
rownames(genetab)<-genetab$Gene
genetab$Gene<-NULL

load_gene <- load_strain %*% t(genetab)

load_gene_df<-as.data.frame(t(load_gene))
info2<-read.csv("genemaf002info.csv",check.names = FALSE)

load_gene_df$Gene<-rownames(load_gene_df)
load_gene_df<-merge(load_gene_df,info2,all = TRUE)
load_gene_df$pos<-load_gene_df$`Genome Fragment`
load_gene_df<-load_gene_df[order(load_gene_df$`Genome Fragment`,load_gene_df$`Order within Fragment`),]

hist(load_gene_df$PC1,xlab = "gene loading",main = "PC1")
plot(load_gene_df$pos,load_gene_df$PC1^2,type = "p")

# outliergenes<-load_gene_df[(load_gene_df$PC1>20)|(load_gene_df$PC1< -20),]
# outliergenes<- outliergenes[order(outliergenes$PC1),]

eggnogdf <- read.table("maf002gene_eggnog.tsv",
                       header = TRUE,
                       check.names = FALSE,
                       fill = TRUE,
                       sep = "\t",
                       quote = "",
                       stringsAsFactors = FALSE)

#anno_outliers<-eggnogdf[eggnogdf$query %in% rownames(outliergenes),]

# library(pracma)
# peaksgene <- findpeaks(load_gene_df$PC1^2, nups = 1, ndowns = 1, minpeakheight = 250)
# peakgene_positions <- peaksgene[, 2]
# load_gene_df$peak<-NA
# load_gene_df$peak[peakgene_positions] <- load_gene_df$PC1[peakgene_positions]

peakgene<-ggplot(load_gene_df,aes(x=pos,y=PC1^2)) +
  geom_point(aes(color = PC1^2 > 300), alpha = 0.5, size = 3) +  # 根据条件设置颜色
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept = 300, color = "blue", linetype = "dashed", linewidth = 1) +  # 添加水平线，y=250
  labs(x = "Accessory Genes (Genome Fragment)", y = "Squared PC1 Loading") +
  theme_classic()+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")+
  scale_y_continuous(breaks = seq(0, 1000, by = 300))
peakgene<-ggplot(load_gene_df, aes(x = PC1^2)) +
  geom_histogram(binwidth = 100,
                 fill = "gray50", 
                 color = "black", 
                 alpha = 0.6) +   # 柱状图
  geom_jitter(aes(y = 0, color = PC1^2 > 300), 
              height = 0.02, 
              size = 3, 
              alpha = 0.5) +    # 点图，颜色按条件分
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_vline(xintercept = 300, 
             color = "blue", 
             linetype = "dashed", 
             linewidth = 1) +  # 阈值线
  labs(x = "Squared PC1 Loading (Accessory Genes)", y = "Count") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"  # 不显示图例，颜色一目了然
  )+
  scale_x_continuous(breaks = seq(0, 1000, by = 300))


anno_outliers2<-eggnogdf[eggnogdf$query %in% load_gene_df[load_gene_df$PC1^2>300, "Gene"],]
colnames(anno_outliers2)[1]<-"Gene"

info2peaks<-info2[info2$Gene %in% load_gene_df[load_gene_df$PC1^2>300, "Gene"],]

genepeakinfo<-merge(anno_outliers2,info2peaks,all = TRUE)
genepeakinfo<-genepeakinfo[order(genepeakinfo$`Order within Fragment`),]
peakinfogen<-genepeakinfo[,c("Genome Fragment","Order within Fragment",
                             "Gene","Annotation","Avg group size nuc",
                             "COG_category","Description")]
peakinfogen$COG_category[is.na(peakinfogen$COG_category) | peakinfogen$COG_category == "-"] <- " "

#write.csv(peakinfogen,"figs/peakgeneinfo.csv",row.names = FALSE,quote = FALSE)
###
eggnogdf <- read.table("maf002gene_eggnog.tsv",
                       header = TRUE,
                       check.names = FALSE,
                       fill = TRUE,
                       sep = "\t",
                       quote = "",
                       stringsAsFactors = FALSE)
anno_gene<-eggnogdf[eggnogdf$query %in% load_gene_df[, "Gene"],]
anno_gene<-anno_gene[,c("query","COG_category")]
library(tidyr)
library(dplyr)
anno_gene <- anno_gene %>%
  mutate(COG_category = strsplit(COG_category, "")) %>%  # 拆成单字母向量
  unnest(COG_category)             
colnames(anno_gene)[1]<-"Gene"
anno_gene <- anno_gene %>% 
  filter(COG_category != "-")

cog_function <- data.frame(
  COG_category = c(
    "J", "A", "K", "L", "B",
    "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O",
    "C", "G", "E", "F", "H", "I", "P", "Q",
    "R", "S"
  ),
  Function = c(
    # INFORMATION STORAGE AND PROCESSING
    "Translation, ribosomal structure and biogenesis",
    "RNA processing and modification",
    "Transcription",
    "Replication, recombination and repair",
    "Chromatin structure and dynamics",
    
    # CELLULAR PROCESSES AND SIGNALING
    "Cell cycle control, cell division, chromosome partitioning",
    "Nuclear structure",
    "Defense mechanisms",
    "Signal transduction mechanisms",
    "Cell wall/membrane/envelope biogenesis",
    "Cell motility",
    "Cytoskeleton",
    "Extracellular structures",
    "Intracellular trafficking, secretion, and vesicular transport",
    "Posttranslational modification, protein turnover, chaperones",
    
    # METABOLISM
    "Energy production and conversion",
    "Carbohydrate transport and metabolism",
    "Amino acid transport and metabolism",
    "Nucleotide transport and metabolism",
    "Coenzyme transport and metabolism",
    "Lipid transport and metabolism",
    "Inorganic ion transport and metabolism",
    "Secondary metabolites biosynthesis, transport and catabolism",
    
    # POORLY CHARACTERIZED
    "General function prediction only",
    "Function unknown"
  ),
  Category = c(
    rep("INFORMATION STORAGE AND PROCESSING", 5),
    rep("CELLULAR PROCESSES AND SIGNALING", 10),
    rep("METABOLISM", 8),
    rep("POORLY CHARACTERIZED", 2)
  )
)

anno_gene <- anno_gene %>%
  left_join(cog_function, by = "COG_category")
genepc1<-load_gene_df[,c("Gene","PC1")]
geneloading_anno<-left_join(anno_gene,genepc1)
geneloading_anno$PC1<-abs(geneloading_anno$PC1)

summary_df <- geneloading_anno %>%
  filter(!is.na(PC1)) %>%
  group_by(COG_category) %>%
  summarise(
    mean_PC1sq = mean(PC1, na.rm = TRUE),
    sd_PC1sq   = sd(PC1, na.rm = TRUE),
    se_PC1sq   = sd(PC1, na.rm = TRUE) / sqrt(n()), 
    n = n(),
    Category = first(Category),  # 或用 unique(Category)[1]
    Function = first(Function),  # 你要是想按 Function 展示位置
    .groups = "drop"
  )

summary_df <- summary_df %>%
  arrange(Category, desc(mean_PC1sq)) %>%
  mutate(Function = factor(Function, levels = unique(Function)))
geneloading_anno <- geneloading_anno %>%
  mutate(Function = factor(Function, levels = levels(summary_df$Function)))

pgene<-ggplot() +
  geom_jitter(
    data = geneloading_anno,
    aes(x = Function, y = PC1, color = Category),
    width = 0.2, alpha = 0.5, size = 2, show.legend = FALSE
  ) +
  geom_errorbar(
    data = summary_df,
    aes(
      x = Function,
      ymin = mean_PC1sq - se_PC1sq,
      ymax = mean_PC1sq + se_PC1sq
    ),
    width = 0.3, linewidth = 0.6
  ) +
  # 平均值横线
  geom_errorbar(
    data = summary_df,
    aes(
      x = Function,
      ymin = mean_PC1sq,
      ymax = mean_PC1sq
    ),
    width = 0.5, linewidth = 1
  )+
  coord_flip() +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_color_npg() +  
  labs(
    x = NULL,
    y = "Absolute Value of PC1 Loadings",
    title = "Accessory Genes"
  ) +
  theme_classic() +
  theme(
    strip.text.y = element_text(size = 12, angle = 0, hjust = 0),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5)
  )+
  scale_x_discrete(limits = rev)+
  geom_hline(yintercept = mean(geneloading_anno$PC1),linewidth = 0.5,linetype = "solid",color="#3B4992FF")
###
my_mean_se <- function(x) {
  m <- mean(x)
  s <- sd(x)
  n <- length(x)
  se <- s / sqrt(n)
  return(c(y = m, ymin = m - se, ymax = m + se))
}
library(ggpubr)
ggplot(geneloading_anno, aes(x = Category, y = PC1, fill = Category)) +
  stat_summary(fun = mean, geom = "bar", width = 0.5, alpha = 0.6, color = "black") +
  stat_summary(fun.data = my_mean_se, geom = "errorbar", width = 0.2, color = "black") +
  geom_jitter(width = 0.1, size = 2, alpha = 0.1, color = "black") +
  stat_compare_means(
    method = "t.test",
    comparisons = list(
      c("CELLULAR PROCESSES AND SIGNALING", "INFORMATION STORAGE AND PROCESSING"),
      c("CELLULAR PROCESSES AND SIGNALING", "METABOLISM"),
      c("INFORMATION STORAGE AND PROCESSING", "METABOLISM"),
      c("INFORMATION STORAGE AND PROCESSING", "POORLY CHARACTERIZED"),
      c("CELLULAR PROCESSES AND SIGNALING", "POORLY CHARACTERIZED"),
      c("METABOLISM", "POORLY CHARACTERIZED")
    ),
    label = "p.format"
  ) +
  theme_light() +
  labs(y = "Projected PC1", x = NULL) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = "none"
  )+scale_fill_npg()
###
genepeaks<-genetab[peakinfogen$Gene,]

library(pheatmap)
pcavec <- read.table("NR369pca.eigenvec")
# remove nuisance column
rownames(pcavec)<-pcavec$V1
pcavec <- pcavec[,-c(1,2)]
colnames(pcavec) <- paste0("PC", 1:ncol(pcavec))
pcaval<-read.table("NR369pca.eigenval")
colnames(pcaval)<-"val"
rownames(pcaval)<-paste0("PC", 1:nrow(pcaval))
pcavec<-as.data.frame(t(t(pcavec)*sqrt(pcaval$val)))

annostrain<-pcavec["PC1"]
colnames(annostrain)<-"strainloading"

annogene<-load_gene_df[load_gene_df$Gene %in% peakinfogen$Gene,c("PC1","Gene")]
annogene$geneloading<-annogene$PC1^2
rownames(annogene)<-annogene$Gene
annogene$Gene<-NULL
annogene$PC1<-NULL

pheatmap(as.matrix(genepeaks),
         show_rownames = FALSE,show_colnames = FALSE,
         annotation_col = annostrain,
         annotation_row = annogene,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE)

compute_ld_r2 <- function(gene1, gene2) {
  n <- length(gene1)  # 样本数
  
  # 计算频率
  P_11 <- sum(gene1 == 1 & gene2 == 1) / n  # 两基因同时出现的频率
  P_1 <- sum(gene1 == 1) / n                # 基因1出现的频率
  P_2 <- sum(gene2 == 1) / n                # 基因2出现的频率
  
  # 避免分母为 0 的情况（边界条件）
  if (P_1 == 0 || P_1 == 1 || P_2 == 0 || P_2 == 1 || (P_1 * (1 - P_1) * P_2 * (1 - P_2)) == 0) {
    return(NA)
  }
  
  # 计算 R^2
  R2 <- (P_11 - P_1 * P_2)^2 / (P_1 * (1 - P_1) * P_2 * (1 - P_2))
  return(R2)
}

# 初始化 LD 矩阵
n_genes <- nrow(genepeaks)  # 基因数量
r2_matrix <- matrix(NA, nrow = n_genes, ncol = n_genes)
rownames(r2_matrix) <- colnames(r2_matrix) <- rownames(genepeaks)

# 计算所有基因对的 LD R^2
for (i in 1:n_genes) {
  for (j in i:n_genes) {  # 只计算上三角部分
    r2_matrix[i, j] <- compute_ld_r2(
      genepeaks[i, ],
      genepeaks[j, ]
    )
    r2_matrix[j, i] <- r2_matrix[i, j]  # 对称矩阵
  }
}

# 检查部分结果
head(r2_matrix)

library(pheatmap)

annodf<-load_gene_df[load_gene_df$Gene %in% peakinfogen$Gene,c("PC1","Gene")]
annodf$PC1<-annodf$PC1^2
rownames(annodf)<-annodf$Gene
annodf$Gene<-NULL

# 绘制热图
pheatmap(
  r2_matrix,
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  main = "LD R^2 Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annodf,
  annotation_col = annodf
)

###all gene LD
# 加载必要的包
library(parallel)

# 初始化基因频率表
genefreq <- as.data.frame(rowSums(genetab) / ncol(genetab))
colnames(genefreq) <- "freq"

# 提取基因总数和基因频率向量
all_genes <- nrow(genetab)  # 基因总数
freq_vector <- genefreq$freq
n <- 369  # 固定样本数

# 初始化 LD 矩阵
allr2_matrix <- matrix(NA, nrow = all_genes, ncol = all_genes)
rownames(allr2_matrix) <- colnames(allr2_matrix) <- rownames(genetab)

# 快速计算函数（带调试信息）
compute_ld_r2 <- function(i, j, freq_vector, genetab, n) {
  if (i == j) {
    return(NA)  # 跳过对角线
  }
  
  gene1 <- genetab[i, ]
  gene2 <- genetab[j, ]
  P_1 <- freq_vector[i]
  P_2 <- freq_vector[j]
  
  # 计算 P_11
  P_11 <- sum(gene1 & gene2) / n
  
  # 计算 R^2
  denominator <- P_1 * (1 - P_1) * P_2 * (1 - P_2)
  if (denominator == 0) {
    return(NA)
  }
  
  R2 <- (P_11 - P_1 * P_2)^2 / denominator
  return(R2)
}

# 并行计算所有基因对的 LD
num_cores <- 15  # 使用固定的核心数
cl <- makeCluster(num_cores)

# 显式导出需要的变量到每个工作节点
clusterExport(cl, list("compute_ld_r2", "freq_vector", "genetab", "n"))

# 并行计算 LD 矩阵
for (i in 1:(all_genes - 1)) {
  # 打印当前进度
  print(paste("Processing gene", i, "out of", all_genes))
  
  # 创建一个索引对列表：对于每个 i，计算所有 j 对
  index_pairs <- cbind(rep(i, all_genes - i), (i + 1):all_genes)
  
  # 导出 index_pairs 给工作节点
  clusterExport(cl, "index_pairs")
  
  # 使用 parSapply 来并行计算每个基因对的 LD 值
  result <- parSapply(cl, 1:nrow(index_pairs), function(idx) {
    compute_ld_r2(index_pairs[idx, 1], index_pairs[idx, 2], freq_vector, genetab, n)
  })
  
  # 将计算结果填充到矩阵的相应位置
  allr2_matrix[i, (i + 1):all_genes] <- result
  allr2_matrix[(i + 1):all_genes, i] <- result
}

# 关闭集群
stopCluster(cl)
##
diag(allr2_matrix) <- 1
###check
checkvalue<-allr2_matrix[genepeakinfo$Gene,genepeakinfo$Gene]

pheatmap(
  checkvalue,
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  main = "LD R^2 Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE
)

#write.csv(allr2_matrix, "figs/allgeneLDr2_matrix.csv", row.names = TRUE)
###prune gene
allr2_matrix<-read.csv("figs/allgeneLDr2_matrix.csv",row.names = 1,
                       check.names = FALSE)
allr2_matrix<-as.matrix(allr2_matrix)
# 设置 LD 阈值
ld_threshold <- 0.8
# 计算 MAF（次要等位基因频率）并存储在 genefreq
maf_values <- genefreq$freq
names(maf_values) <- rownames(genefreq)  # 将基因名称作为名字
# 初始化保留基因列表
selected_genes <- colnames(allr2_matrix)
# 遍历 LD 矩阵的上三角部分，筛选出高 LD 基因对
for (i in 1:(ncol(allr2_matrix) - 1)) {
  for (j in (i + 1):ncol(allr2_matrix)) {
    # 检查 LD 值是否超过阈值
    if (!is.na(allr2_matrix[i, j]) && allr2_matrix[i, j] > ld_threshold) {
      gene_i <- colnames(allr2_matrix)[i]
      gene_j <- colnames(allr2_matrix)[j]
      
      # 比较 MAF 值，移除 MAF 较低的基因
      if (maf_values[gene_i] < maf_values[gene_j]) {
        selected_genes <- setdiff(selected_genes, gene_i)
      } else {
        selected_genes <- setdiff(selected_genes, gene_j)
      }
    }
  }
}
pruned_genetable <- genetab[selected_genes, ]
pruned_R2<- allr2_matrix[selected_genes,selected_genes]
##genepca
library(dplyr)
vcf_fixed <- data.frame(
  CHROM = "1",                      # 假设所有基因属于同一染色体
  POS = 1:nrow(pruned_genetable),      # 使用基因的行号作为位置
  ID = rownames(pruned_genetable),     # 基因名称作为位点 ID
  REF = "A",                        # 参考等位基因（默认 A）
  ALT = "T",                        # 替代等位基因（默认 T）
  QUAL = ".",                       # 缺失值（可选）
  FILTER = "PASS",                  # 标记为通过过滤（默认 PASS）
  INFO = ".",                       # 信息列（可空）
  FORMAT = "GT"                     # 基因型格式
)
vcf_data <- cbind(vcf_fixed, pruned_genetable)
# write.table(
#   vcf_data,
#   "figs/pruned_genetable.vcf",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )
genepca<-read.table("figs/genepca.eigenvec")
geneval<-read.table("figs/genepca.eigenval")
rownames(genepca)<-genepca$V1
genepca <- genepca[,-c(1,2)]
colnames(genepca) <- paste0("PC", 1:ncol(genepca))

genepca<-as.data.frame(t(t(genepca)*sqrt(geneval$V1)))
plot(genepca$PC3,pca$PC1,type = "p",
     xlab = "PC3 based on genes",ylab = "PCA based on SNPs")

###enrichment
library(tidytree)
library(dplyr)
library(stringr)
library(clusterProfiler)
l<-unique(anno_outliers$query)
data<-eggnogdf[,c("query","GOs")]
data_cog <- data[data$GOs!="-",]
data_cog <- data_cog %>%
  mutate(GOs = str_split(GOs,",")) %>%
  unnest(GOs)

data2gene_list <- split(data_cog$query, data_cog$GOs)
cog2gene <- stack(data2gene_list)
colnames(cog2gene) <- c("GENE", "TERM")
cog2gene<-cog2gene[,c("TERM","GENE")]

ego <- enricher(
  l,
  TERM2GENE = cog2gene
)
barplot(ego)
dotplot(ego)
###meanvalues for each cog
load_gene_df$query<-rownames(load_gene_df)
cogdf<-merge(load_gene_df,eggnogdf,all = TRUE)
cogdf<-na.omit(cogdf[,c("PC1","COG_category")])
cogdf<-cogdf[cogdf$COG_category!="-",]
cogdf <- cogdf %>%
  mutate(COG_category = str_split(COG_category,"")) %>%
  unnest(COG_category)
cogdf$PC1<-abs(cogdf$PC1)

cogsum <- cogdf %>%
  group_by(COG_category) %>%
  summarise(mean = mean(PC1),
            count = n(),
            se = sd(PC1)/sqrt(n()))

cogsum<-cogsum[order(cogsum$mean),]
cogsum$COG_category<-factor(cogsum$COG_category,levels = unique(cogsum$COG_category))

library(ggplot2)
ggplot(cogsum,aes(x=COG_category))+
  geom_bar(aes(y=mean),stat = "identity",fill="lightblue",alpha=0.7)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.2)+
  geom_text(aes(y=mean,label = count),vjust= -7,color="black")+
  labs(title = "COG category Mean, Count, and Standard Error",y = "Absolute Mean Loading")+
  theme_minimal()+
  geom_hline(yintercept = mean(cogsum$mean),linetype="dashed",color="red",linewidth=1)

ggplot(cogdf,aes(x=COG_category,y=PC1))+
  geom_violin(fill="lightblue",alpha=0.7)+
  geom_boxplot(width=0.1,fill="white")+
  theme_minimal()

###snp
geno<-read.table("NR369maf002.matrix_01",check.names = FALSE)
load_snp <- load_strain %*% t(geno)
load_snp_squ<-load_snp^2

load_snp_df<-as.data.frame(t(load_snp))
load_snp_df$POS<-as.numeric(rownames(load_snp_df))

###
geneloading_anno$Group <- "Accessory Genes"
meangene$Group <- "Core Genes"
meangene$PC1<-meangene$avg_PC1
g<-geneloading_anno[,c("PC1","Group")]
s<-meangene[,c("PC1","Group")]
load_snp_df$Group <- "Core SNPs"
cs<-load_snp_df[,c("PC1","Group")]
cs$PC1<-abs(cs$PC1)
df_all <- bind_rows(g, s,cs)

pgenesnp<-ggplot(df_all, aes(x = PC1, fill = Group, color=Group)) +
  geom_density(alpha = 0.5) +
  theme_light() +
  scale_fill_aaas()+
  scale_color_aaas()+
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1, "lines"),
    legend.position = c(0.6,0.8)
  ) + ylab("Density") + xlab("Absolute PC1 Loading")+
  geom_vline(xintercept = mean(geneloading_anno$PC1),color="#3B4992FF", linetype = "dashed")+
  geom_vline(xintercept = mean(meangene$PC1),color="#EE0000FF", linetype = "dashed")+
  geom_vline(xintercept = mean(cs$PC1),color="#008B45FF", linetype = "dashed")

###
squdf<-as.data.frame(t(load_snp_squ))
squdf$POS<-as.numeric(rownames(squdf))
library(ggplot2)
p1<-ggplot(squdf,aes(x=POS,y=PC1,color=POS))+
  geom_point()+
  theme_classic()+
  xlab("Position (Mb)")+
  ylab("Squared PC1 Loading")+
  scale_x_continuous(labels = c("0","1","2","3","4","5"))+
  scale_color_gradientn(colours = c("firebrick", "springgreen4", 
                                    "blue3", "turquoise3","darkorchid2", "gold2","orange"),
                        guide="none")+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16))
library(cowplot)
plot_grid(p1,p2,p3,p4,p5,align = "hv",ncol = 1)

supL<-ggplot(load_snp_df,aes(x=POS,y=PC1))+
  geom_point()+
  theme_classic()+
  xlab("Position (Mb)")+
  ylab("PC1 Loading")+
  scale_x_continuous(labels = c("0","1","2","3","4","5"))+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16))
###
# library(pracma)
# peaks <- findpeaks(squdf$PC1, nups = 1, ndowns = 1, minpeakheight = 250)
# 
# plot(squdf$POS, squdf$PC1, type = "p", col = "black", 
#      main = "Peak Detection on Squared PC1 Loading",
#      xlab = "Position (Mb)", ylab = "Squared PC1 Loading")
# points(squdf$POS[peaks[,2]], squdf$PC1[peaks[,2]], col = "red", pch = 19)  # 标记红色的峰值点
# 
# peak_positions <- peaks[, 2]
# squdf$peak<-NA
# squdf$peak[peak_positions] <- squdf$PC1[peak_positions]
peaksnp<-ggplot(squdf,aes(x=POS,y=PC1)) +
  geom_point(aes(color = PC1 > 300), alpha = 0.5, size = 3) +  # 根据条件设置颜色
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept = 300, color = "blue", linetype = "dashed", linewidth = 1) +  # 添加水平线，y=250
  labs(x = "Core SNPs Position (Mb)", y = "Squared PC1 Loading") +
  scale_x_continuous(labels = c("0","1","2","3","4","5"))+
  theme_classic()+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")+
  geom_vline(xintercept = 300000)+geom_vline(xintercept = 400000)+
  annotate("text", x = 200000, y = 750, label = "1", size = 6,fontface = "bold")+
  geom_vline(xintercept = 1520000)+geom_vline(xintercept = 1620000)+
  annotate("text", x = 1420000, y = 750, label = "2", size = 6,fontface = "bold")+
  geom_vline(xintercept = 1880000)+geom_vline(xintercept = 1980000)+
  annotate("text", x = 2080000, y = 750, label = "3", size = 6,fontface = "bold")+
  geom_vline(xintercept = 5080000)+geom_vline(xintercept = 5180000)+
  annotate("text", x = 4980000, y = 750, label = "4", size = 6,fontface = "bold")+
  scale_y_continuous(breaks = seq(0, 1000, by = 300))

library(cowplot)
fig1c<-plot_grid(peaksnp,peakgene,
          align = "hv",nrow = 1,
          rel_widths = c(2, 1))
ggsave("/Users/47liu/Documents/paper/peer review/figs/figure1c.pdf",
       plot = fig1c,
       width = 16, height = 5,
       units = "in",
       dpi = 600,
       bg = "transparent",
       useDingbats = FALSE  
)
ggsave("/Users/47liu/Documents/paper/peer review/figs/fig1c.png", fig1c,
       width = 16, height = 5, units = "in", dpi = 600)
###tree1 vs tree2
ggplot(squdf, aes(x = POS, y = PC1)) +
  geom_point(data = subset(squdf, PC1 > 300), aes(color = "above_300"), alpha = 0.5, size = 3) +
  geom_point(data = subset(squdf, PC1 < 10), aes(color = "below_100"), alpha = 0.5, size = 3) +
  geom_point(data = subset(squdf, PC1 >= 10 & PC1 <= 300), aes(color = "default"), alpha = 0.5, size = 3) +
  
  # 自定义颜色
  scale_color_manual(values = c(
    "above_300" = "red",  # PC1 > 300 的点为红色
    "below_100" = "blue", # PC1 < 100 的点为蓝色
    "default" = "black"   # 其他情况为黑色
  )) +
  
  # 添加水平线
  geom_hline(yintercept = 300, color = "red", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 10, color = "blue", linetype = "dashed", linewidth = 1) +
  
  # 设置标签和主题
  labs(x = "SNP Position (Mb)", y = "Squared PC1 Loading") +
  scale_x_continuous(labels = c("0", "1", "2", "3", "4", "5")) +
  scale_y_continuous(breaks = seq(0, 1000, by = 300)) +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 1),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 20),
    legend.position = "none"
  )
tree1pos<-squdf[squdf$PC1>300,"POS"]
tree2pos<-squdf[squdf$PC1<10,"POS"]
##
#writeLines(as.character(tree1pos), con = "tree1pos.txt")
#writeLines(as.character(tree2pos), con = "tree2pos.txt")
library(ggtree)
library(treeio)
tree1<-read.newick("tree1pos.fasttree")
tree2<-read.newick("tree2pos.fasttree")

pca <- read.table("NR369pca.eigenvec")
rownames(pca)<-pca$V1
pca <- pca[,-c(1,2)]
colnames(pca) <- paste0("PC", 1:ncol(pca))
val<-read.table("NR369pca.eigenval")
colnames(val)<-"val"
rownames(val)<-paste0("PC", 1:nrow(val))
pca<-as.data.frame(t(t(pca)*sqrt(val$val)))

pca$label<-rownames(pca)
tree1.a<-full_join(tree1,pca,by="label")
tree2.a<-full_join(tree2,pca,by="label")
library(viridis)
ggtree(tree1.a,layout = "ape")+
  geom_tippoint(aes(color = PC1)) +  
  scale_color_viridis(
    option = "H",  # 使用 Viridis 渐变色
    name = "PC1",  # 图例标题
    breaks = seq(-0.3, 0.4, by = 0.2),  # 设置图例断点
    guide = guide_colorbar(reverse = TRUE)
  )+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(2, "lines"),
        legend.position = "right")


ggtree(tree2.a,layout = "ape")+
  geom_tippoint(aes(color = PC1)) +  
  scale_color_viridis(
    option = "H",  # 使用 Viridis 渐变色
    name = "PC1",  # 图例标题
    breaks = seq(-0.3, 0.4, by = 0.2),  # 设置图例断点
    guide = guide_colorbar(reverse = TRUE)
  )+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(2, "lines"),
        legend.position = "right")


###ld snp
snplist<-squdf[squdf$PC1>300,"POS"]
snppeak<-geno[as.character(snplist),]

compute_ld_r2 <- function(gene1, gene2) {
  n <- length(gene1)  # 样本数
  
  # 计算频率
  P_11 <- sum(gene1 == 1 & gene2 == 1) / n  # 两基因同时出现的频率
  P_1 <- sum(gene1 == 1) / n                # 基因1出现的频率
  P_2 <- sum(gene2 == 1) / n                # 基因2出现的频率
  
  # 避免分母为 0 的情况（边界条件）
  if (P_1 == 0 || P_1 == 1 || P_2 == 0 || P_2 == 1 || (P_1 * (1 - P_1) * P_2 * (1 - P_2)) == 0) {
    return(NA)
  }
  
  # 计算 R^2
  R2 <- (P_11 - P_1 * P_2)^2 / (P_1 * (1 - P_1) * P_2 * (1 - P_2))
  return(R2)
}

# 初始化 LD 矩阵
n_snp <- nrow(snppeak)  # snp数量
snpr2_matrix <- matrix(NA, nrow = n_snp, ncol = n_snp)
rownames(snpr2_matrix) <- colnames(snpr2_matrix) <- rownames(snppeak)

# 计算所有基因对的 LD R^2
for (i in 1:n_snp) {
  print(i)
  for (j in i:n_snp) {  # 只计算上三角部分
    snpr2_matrix[i, j] <- compute_ld_r2(
      snppeak[i, ],
      snppeak[j, ]
    )
    snpr2_matrix[j, i] <- snpr2_matrix[i, j]  # 对称矩阵
  }
}

library(pheatmap)

annodf2<-squdf[squdf$POS %in% snplist,c("PC1","POS")]
rownames(annodf2)<-annodf2$POS

# 绘制热图
pheatmap(
  snpr2_matrix,
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  main = "LD R^2 Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annodf2,
  annotation_col = annodf2
)

###ld snp and gene
peakdf<-rbind(genepeaks,snppeak)

# 初始化 LD 矩阵
n_variant <- nrow(peakdf)  
Vr2_matrix <- matrix(NA, nrow = n_variant, ncol = n_variant)
rownames(Vr2_matrix) <- colnames(Vr2_matrix) <- rownames(peakdf)

# 计算所有基因对的 LD R^2
for (i in 1:n_variant) {
  print(i)
  for (j in i:n_variant) {  # 只计算上三角部分
    Vr2_matrix[i, j] <- compute_ld_r2(
      peakdf[i, ],
      peakdf[j, ]
    )
    Vr2_matrix[j, i] <- Vr2_matrix[i, j]  # 对称矩阵
  }
}
annodf$type<-"Gene"
annodf2$type<-"SNP"
annodf2$POS<-NULL

annoall<-rbind(annodf,annodf2)
annoall$position<-as.numeric(rownames(annoall))

pheatmap(
  Vr2_matrix,
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  main = "LD R^2 Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annoall,
  annotation_col = annoall,
  legend = FALSE,
  annotation_legend = FALSE
)
#write.csv(Vr2_matrix, "figs/Vr2_matrix.csv", row.names = TRUE)
# annocoldf<-pca["PC1"]
# pheatmap(
#   as.matrix(peakdf),
#   cluster_rows = TRUE,  
#   cluster_cols = TRUE,  
#   show_rownames = FALSE,
#   show_colnames = FALSE,
#   annotation_col = annocoldf,
#   annotation_row = annoall
# )
###add more information
Vr2_matrix<-as.matrix(read.csv("figs/Vr2_matrix.csv",row.names = 1,check.names = FALSE))
###remove group_205
Vr2_matrix <- Vr2_matrix[!rownames(Vr2_matrix) %in% "group_205", 
                         !colnames(Vr2_matrix) %in% "group_205"]
annoall<-annoall[!rownames(annoall) %in% "group_205",]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggsci)
# 确保 position 列中没有 NA 的最小值和最大值被正确计算
position_min <- 0
position_max <- 5248520
# 定义符合正刊的配色方案
nature_colormap <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
# 分类颜色 (type)
type_colors <- c("Gene" = pal_cosmic("hallmarks_light")(2)[1], "SNP" = pal_cosmic("hallmarks_light")(2)[2])
# PC1 颜色柔和
pc1_colors <- colorRamp2(
  c(min(annoall$PC1, na.rm = TRUE), 
    mean(annoall$PC1, na.rm = TRUE), 
    max(annoall$PC1, na.rm = TRUE)), 
  c("#E5F5F9", "#99D8C9", "#2CA25F") # 浅蓝到深绿
)
# Position 颜色
position_colors <- colorRamp2(
  c(position_min, position_min + (position_max - position_min) / 3, 
    position_min + 2 * (position_max - position_min) / 3, position_max), # 4 个色阶
  c("#FEE6CE", "#FDAE6B", "#E6550D", "#A63603") # 从浅橙到深棕的四个颜色
)
# 构建列注释（带图例）
ha <- rowAnnotation(
  df = annoall, 
  col = list(
    type = type_colors,
    PC1 = pc1_colors,
    position = position_colors
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    type = list(title = "Type",
                title_gp = gpar(fontsize = 14),  
                labels_gp = gpar(fontsize = 14)),     
    PC1 = list(title = expression(PC1^2),
               at=seq(300,1000,by=300),
               title_gp = gpar(fontsize = 14),  
               labels_gp = gpar(fontsize = 14)), 
    position = list(
      title = "Position (Mb)",  # 修改 position 图例的标题
      at = seq(position_min, position_max, by = 1000000),  # 设置 position 的显示范围和间隔
      labels = c("0", "1", "2", "3", "4", "5"),  # 设置图例标签
      title_gp = gpar(fontsize = 14),  
      labels_gp = gpar(fontsize = 14)) 
  )
)
# 绘制热图
Heatmap(
  Vr2_matrix,
  name = "LD R^2", # Heatmap title
  col = nature_colormap, # Use Nature-style color mapping
  cluster_rows = TRUE, # Perform row clustering
  cluster_columns = TRUE, # Perform column clustering
  show_row_names = FALSE, # Do not show row names
  show_column_names = FALSE, # Do not show column names
  left_annotation = ha, # Add annotation
  show_row_dend = FALSE,
  column_dend_height = unit(40, "mm")
)

###cut heatmap
# 创建热图对象，但不立即绘制
ht <- Heatmap(
  Vr2_matrix,
  name = "LD R^2", # Heatmap title
  col = nature_colormap, # Use Nature-style color mapping
  cluster_rows = TRUE, # Perform row clustering
  cluster_columns = TRUE, # Perform column clustering
  show_row_names = FALSE, # Do not show row names
  show_column_names = FALSE, # Do not show column names
  left_annotation = ha, # Add annotation
  show_row_dend = FALSE,
  column_dend_height = unit(40, "mm")
)

# 绘制热图，并生成热图对象
ht_draw <- draw(ht)

# 提取行和列的顺序
row_order <- row_order(ht_draw) # 行顺序
col_order <- column_order(ht_draw) # 列顺序

# 按顺序重新排列矩阵
Vr2_matrix_ordered <- Vr2_matrix[row_order, col_order]

# 提取行的聚类树
row_dend <- row_dend(ht_draw)  # 提取行的 dendrogram

# 提取列的聚类树
column_dend <- column_dend(ht_draw)  # 提取列的 dendrogram

# 提取列的聚类树
column_dend <- column_dend(ht_draw)

# 将列聚类树切分为 7 组
column_groups <- cutree(as.hclust(column_dend), k = 7)

# 使用 dendextend 美化聚类树
library(dendextend)
colored_column_dend <- color_branches(column_dend, k = 7)  # 对分组进行着色

# 使用 rainbow() 为分组设置颜色
group_colors <- rainbow(7)  # 生成 7 种颜色
names(group_colors) <- as.character(1:7)  # 将颜色与分组标签对应

# 使用 dendextend 美化聚类树
colored_column_dend <- column_dend %>% 
  color_branches(k = 7, col = group_colors)  # 应用分组颜色

# 绘制美化后的聚类树
plot(colored_column_dend, main = "Column Dendrogram with 7 Groups (Rainbow Colors)")

column_info <- data.frame(
  Column = colnames(Vr2_matrix),  # 列名称
  Group = column_groups,          # 分组编号
  Color = group_colors[as.character(column_groups)]  # 对应的颜色
)

# 定义映射关系
group_mapping <- c("5" = "1", "1" = "Other 40 blocks", "2" = "3", 
                   "3" = "Other 40 blocks", "4" = "Other 40 blocks", "6" = "2", "7" = "4")

# 映射 Group 列到 Mapped_Group
column_info$Mapped_Group <- group_mapping[as.character(column_info$Group)]

###3a,3b,3c
a <-rownames(rightdf)[1:26]
b<-rownames(rightdf)[27:37]
c<-rownames(rightdf)[38:nrow(rightdf)]
column_info[a,"Mapped_Group"]<-"3a"
column_info[b,"Mapped_Group"]<-"3b"
column_info[c,"Mapped_Group"]<-"3c"

library(viridis)
# 定义颜色映射规则（为 Mapped_Group 分配颜色）
mapped_group_colors <- viridis(length(unique(column_info$Mapped_Group)), option = "D") # 冷暖色调
names(mapped_group_colors) <- unique(column_info$Mapped_Group)

# 创建 HeatmapAnnotation 用于顶部颜色条
mapped_group_bar <- HeatmapAnnotation(
  Mapped_Group = factor(column_info$Mapped_Group, levels = sort(unique(column_info$Mapped_Group))), # 确保按顺序排列
  col = list(Mapped_Group = mapped_group_colors), # 为分组定义颜色
  annotation_legend_param = list(title = "Block",
                                 title_gp = gpar(fontsize = 14),  
                                 labels_gp = gpar(fontsize = 14)), # 移除图例标题
  show_annotation_name = FALSE # 隐藏注释名称
)

#pdf("/Users/47liu/Documents/paper/peer review/figs/figure1d.pdf", width = 16, height = 6)
png("/Users/47liu/Documents/paper/peer review/figs/fig1d.png",
    width = 16, height = 6, units = "in", res = 600)
# 绘制热图并添加顶部颜色条
Heatmap(
  Vr2_matrix,
  name = "LD R^2", # Heatmap title
  col = nature_colormap, # Use Nature-style color mapping
  cluster_rows = TRUE, # Perform row clustering
  cluster_columns = TRUE, # Perform column clustering
  show_row_names = FALSE, # Do not show row names
  show_column_names = FALSE, # Do not show column names
  #bottom_annotation = mapped_group_bar,  # 在热图顶部添加颜色条
  left_annotation = ha, # Add annotation
  show_row_dend = FALSE,
  column_dend_height = unit(20, "mm"),
  heatmap_legend_param = list(
    title = expression("LD " * R^2),  # 设置图例标题
    title_gp = gpar(fontsize = 14),  
    labels_gp = gpar(fontsize = 14) 
  )
)
dev.off()

###
#annotation
anno1<-read.csv("NR369maf002anno.txt",sep = "|",header = TRUE,fill = TRUE)
anno2<-read.csv("MM_yezscq0m.emapper.annotations.tsv",sep = "\t")
library(dplyr)
snp_info <-left_join(anno1, anno2)
colnames(snp_info)[1]<-"position"
allsnpinfo<-left_join(annoall,snp_info)
allsnpinfo<-allsnpinfo[c(52:nrow(allsnpinfo)),]
colnames(allsnpinfo)[3]<-"label"
allsnpinfo$label<-as.character(allsnpinfo$label)

info2<-read.csv("genemaf002info.csv",check.names = FALSE)
eggnogdf <- read.table("maf002gene_eggnog.tsv",
                       header = TRUE,
                       check.names = FALSE,
                       fill = TRUE,
                       sep = "\t",
                       quote = "",
                       stringsAsFactors = FALSE)
colnames(eggnogdf)[1]<-"Gene"
gene_info<-left_join(info2,eggnogdf)
genehits<-gene_info[gene_info$Gene %in% rownames(annoall),]
colnames(genehits)[1]<-"label"

allinfo<-annoall
allinfo$label<-rownames(annoall)

step1<-merge(allinfo,allsnpinfo)
colnames(step1)[8]<-"SNP_Annotation"

step2<-merge(allinfo,genehits)
colnames(step2)[6]<-"Gene_Annotation"

hitsinfo<-merge(step1, step2,all = TRUE)
hitsinfokeep<-hitsinfo[,c("label","PC1","type","position",
                          "SNP_Annotation","Gene_Annotation","Gene_ID","product",
                          "COG_category","Description")]

column_info$label<-column_info$Column
hitsinfokeep_sortedblocks<-merge(hitsinfokeep,column_info,all = TRUE)
hitsinfokeep_sortedblocks <- hitsinfokeep_sortedblocks[match(rownames(Vr2_matrix_ordered), hitsinfokeep_sortedblocks$label), ]
hitsinfokeep_sortedblocks<-hitsinfokeep_sortedblocks[,c("label","PC1","type","position","SNP_Annotation",
                                                        "Gene_Annotation","Gene_ID","product","COG_category","Description",
                                                        "Mapped_Group")]
colnames(hitsinfokeep_sortedblocks)[11]<-"Block"
rownames(hitsinfokeep_sortedblocks)<-hitsinfokeep_sortedblocks$label
library(writexl)

hitsinfokeep_sortedblocks$Gene_ID <- gsub("^gene-", "", hitsinfokeep_sortedblocks$Gene_ID)
hitsinfokeep_sortedblocks$Gene_ID <- replace(hitsinfokeep_sortedblocks$Gene_ID, 
                                             hitsinfokeep_sortedblocks$Gene_ID == "null", 
                                             NA)
hitsinfokeep_sortedblocks$COG_category <- replace(hitsinfokeep_sortedblocks$COG_category, 
                                             hitsinfokeep_sortedblocks$COG_category == "-", 
                                             NA)
library(readxl)
hitsanno<-read_excel("figs/20250107heatmapAnnotation.xlsx")
hitsanno <- hitsanno[hitsanno$label != "group_205", ]
hitsanno_sorted <- hitsanno[match(hitsinfokeep_sortedblocks[[1]], hitsanno[[1]]), ]

#rownames(hitsinfokeep_sortedblocks)<-c(1:nrow(hitsinfokeep_sortedblocks))

hitsinfokeep_sortedblocks[252:nrow(hitsinfokeep_sortedblocks),"other40blocks"]<-factor(c(rep(1,9),rep(2,5),3,4,4,rep(5,4),rep(6,5),rep(7,7),c(8:13),14,14,15,rep(16,6),
                                                                                         17,17,18,18,18,19,19,19,rep(20,7),21,21,22,23,24,25,rep(26,4),27,27,c(28:40)))

#write_xlsx(hitsinfokeep_sortedblocks, path = "figs/heatmapAnnotation.xlsx")
###look blocks
block_assignment <- column_info$Mapped_Group # 提取 Block 信息
Vr2_matrix<-Vr2_matrix[column_info$Column,column_info$Column]
unique_blocks <- unique(block_assignment) # 获取唯一的 Block 值
# 定义主图的全局值范围（根据主图计算或手动指定）
global_min <- min(Vr2_matrix, na.rm = TRUE)
global_max <- max(Vr2_matrix, na.rm = TRUE)
# 定义配色方案（与主热图一致）
library(RColorBrewer)
block_colormap <- colorRamp2(
  seq(global_min, global_max, length.out = 100), # 强制颜色范围固定为主热图的范围
  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
)
#1234
# 提取 Block 1234 的行列索引
block_1_idx <- which(grepl("4", block_assignment))
block_1_idx <- which(block_assignment=="4")
# 提取 Block 1234 的子矩阵
block_1_matrix <- Vr2_matrix[block_1_idx, block_1_idx]
block_1_matrix <- block_1_matrix[hitsinfokeep_sortedblocks[grepl("4", hitsinfokeep_sortedblocks$Block),"label"],
                                 hitsinfokeep_sortedblocks[grepl("4", hitsinfokeep_sortedblocks$Block),"label"]]
block_1_matrix <- block_1_matrix[hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="4","label"],
                                 hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="4","label"]]
block_1_annotation <- annoall[rownames(block_1_matrix), ] # 根据行索引提取注释数据
# 定义左侧注释条（rowAnnotation）
left_annotation <- rowAnnotation(
  df = block_1_annotation, # 注释数据框
  col = list(
    type = type_colors,       # 定义 type 注释的颜色映射
    PC1 = pc1_colors,         # 定义 PC1 注释的颜色映射
    position = position_colors # 定义 position 注释的颜色映射
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    type = list(title = "Type",
                title_gp = gpar(fontsize = 13),  
                labels_gp = gpar(fontsize = 12)),     
    PC1 = list(title = expression(PC1^2),
               at=seq(300,1000,by=300),
               title_gp = gpar(fontsize = 13),  
               labels_gp = gpar(fontsize = 12)), 
    position = list(
      title = "Position (Mb)",  # 修改 position 图例的标题
      at = seq(position_min, position_max, by = 1000000),  # 设置 position 的显示范围和间隔
      labels = c("0", "1", "2", "3", "4", "5"),  # 设置图例标签
      title_gp = gpar(fontsize = 13),  
      labels_gp = gpar(fontsize = 12)) 
  )
)
rightdf<-hitsinfokeep_sortedblocks[grepl("4", hitsinfokeep_sortedblocks$Block),
                                   c("Gene_ID","SNP_Annotation","COG_category")]
rightdf<-hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="4",
                                   c("Gene_ID","SNP_Annotation","COG_category")]
r<-rownames(rightdf)
b3<-hitsanno_sorted[grepl("4",hitsanno_sorted$Block), ]
b3<-hitsanno_sorted[hitsanno_sorted$Block=="4", ]
library(dplyr)
b3 <- b3 %>%
  mutate(
    `old locus tag(NC_012731.1)` = case_when(
      Gene_ID == "KP1_RS30365" ~ Gene_ID,
      Gene_ID == "KP1_RS30360" ~ Gene_ID,
      TRUE ~ `old locus tag(NC_012731.1)`  # 保持原值
    )
  )
rightdf$Gene_ID<-b3$`old locus tag(NC_012731.1)`
rownames(rightdf)<-r
rightdf <- rightdf[rownames(block_1_matrix), ]
# 确保字符型列是因子（因子便于颜色映射）
rightdf$Gene_ID <- factor(rightdf$Gene_ID)
rightdf$SNP_Annotation <- factor(rightdf$SNP_Annotation)
rightdf$COG_category <- factor(rightdf$COG_category)
# 定义颜色映射规则（基于因子级别）
gene_id_colors <- setNames(
  pal_lancet("lanonc")(length(levels(rightdf$Gene_ID))), # 提取对应数量的颜色
  levels(rightdf$Gene_ID)
)
# gene_id_colors <- setNames(
#   rainbow(27),  # 使用rainbow生成27种颜色
#   levels(rightdf$Gene_ID)
# )
snp_annotation_colors <- setNames(
  pal_jama("default")(length(levels(rightdf$SNP_Annotation))), # 提取对应数量的颜色
  levels(rightdf$SNP_Annotation)
)
cog_category_colors <- setNames(
  pal_npg("nrc")(length(levels(rightdf$COG_category))), # 提取对应数量的颜色
  levels(rightdf$COG_category)
)
# 创建右侧注释条
right_annotation <- rowAnnotation(
  df = rightdf, # 注释数据框
  col = list(
    Gene_ID = gene_id_colors,          # Gene_ID 的颜色映射
    SNP_Annotation = snp_annotation_colors,    # SNP_Annotation 的颜色映射
    COG_category = cog_category_colors  # COG_category 的颜色映射
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Gene_ID = list(title = "Gene",
                title_gp = gpar(fontsize = 13),  
                labels_gp = gpar(fontsize = 12),
                ncol=1),     
    SNP_Annotation = list(title = "SNP",
               title_gp = gpar(fontsize = 13),  
               labels_gp = gpar(fontsize = 12)), 
    COG_category = list(
      title = "COG",  
      title_gp = gpar(fontsize = 13),  
      labels_gp = gpar(fontsize = 12),
      ncol=1) 
  )
)
ht1 <- Heatmap(
  block_1_matrix,
  name = "LD R^2", 
  col = block_colormap, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  column_dend_height = unit(20, "mm"),
  show_row_dend = FALSE,
  left_annotation = left_annotation,
  right_annotation = right_annotation,
  heatmap_legend_param = list(
    title = expression("LD " * R^2),  # 设置图例标题
    title_gp = gpar(fontsize = 13),  
    labels_gp = gpar(fontsize = 12) 
  )
)
draw(
  ht1,
  heatmap_legend_side = "bottom",          
  annotation_legend_side = "bottom",
  merge_legend = TRUE,
  padding = unit(c(0.2, 0, 1, 0), "cm") 
)
grid.text(
  "Block 4",                     # 替换为您的标题
  x = unit(0.5, "npc"),                   # 标题水平居中
  y = unit(0.97, "npc"),                  # 标题在顶部，稍微留出间距
  gp = gpar(fontsize = 16, fontface = "bold") # 字体大小和加粗
)
##
#pcvalues<-pca[order(pca$PC1),"ID"]
block_1_0_1matrix<-as.matrix(peakdf[rownames(block_1_matrix),])
#row_clustering <- hclust(dist(block_1_0_1matrix))
col_clustering <- hclust(dist(t(block_1_0_1matrix))) 
block_1_0_1matrix <- apply(block_1_0_1matrix, c(1, 2), function(x) factor(x, levels = c(0, 1)))

#strain_colors <- c("#FF0000", "#800080", "#8A2BE2", "#7FBF7F", "#00FF00")
variant_colors <- c("0" = "#4575B4", "1" = "#D73027")

# 优化颜色映射
annostrain <- annostrain[colnames(block_1_0_1matrix), , drop = FALSE]
color_mapping <- circlize::colorRamp2(
  seq(min(annostrain$strainloading), max(annostrain$strainloading), length.out = 369),  # 创建渐变范围
  viridis(369, option = "H")  # 使用 viridis 的 "H" 配色
)
# color_mapping <- colorRamp2(
#   seq(min(annostrain$strainloading), 
#       max(annostrain$strainloading), 
#       length.out = length(strain_colors)), 
#   strain_colors
# )

strains_annotation <- HeatmapAnnotation(
  Strain = annostrain$strainloading,  
  col = list(Strain = color_mapping),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Strain = list(title = "Strain",
                  at=seq(-0.3,0.4,by=0.2),
                title_gp = gpar(fontsize = 13),  
                labels_gp = gpar(fontsize = 12))
  )
)

ht1.1 <- Heatmap(
  block_1_0_1matrix,
  name = "Variants", 
  cluster_rows = FALSE,  
  cluster_columns = as.dendrogram(col_clustering),
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  #show_row_dend = FALSE,
  left_annotation = left_annotation,
  right_annotation = right_annotation,
  top_annotation = strains_annotation,
  col = variant_colors,  # 应用颜色映射
  heatmap_legend_param = list(
    title = "Variant",  # 设置图例标题
    title_gp = gpar(fontsize = 13),  
    labels_gp = gpar(fontsize = 12) 
  )
)

draw(
  ht1.1,
  heatmap_legend_side = "bottom",          
  annotation_legend_side = "bottom",
  merge_legend = TRUE,
  padding = unit(c(0.2, 0, 1, 0), "cm") 
)
# grid.text(
#   "Block 3",                     # 替换为您的标题
#   x = unit(0.5, "npc"),                   # 标题水平居中
#   y = unit(0.97, "npc"),                  # 标题在顶部，稍微留出间距
#   gp = gpar(fontsize = 16, fontface = "bold") # 字体大小和加粗
# )
#####
# #2
# block_2_idx <- which(block_assignment == 2)
# block_2_matrix <- Vr2_matrix[block_2_idx, block_2_idx]
# block_2_matrix <- block_2_matrix[hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="2","label"],
#                                  hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="2","label"]]
# block_2_annotation <- annoall[rownames(block_2_matrix), ] # 根据行索引提取注释数据
# # 定义左侧注释条（rowAnnotation）
# left_annotation <- rowAnnotation(
#   df = block_2_annotation, # 注释数据框
#   col = list(
#     type = type_colors,       # 定义 type 注释的颜色映射
#     PC1 = pc1_colors,         # 定义 PC1 注释的颜色映射
#     position = position_colors # 定义 position 注释的颜色映射
#   ),
#   show_annotation_name = FALSE,
#   annotation_legend_param = list(
#     type = list(title = "Type"),     # 修改 type 图例的标题
#     PC1 = list(title = "PC1^2"), # 修改 PC1 图例的标题
#     position = list(title = "Position") # 修改 position 图例的标题
#   )
# )
# rightdf<-hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block==2,
#                                    c("Gene_ID","SNP_Annotation","COG_category")]
# rightdf <- rightdf[rownames(block_2_matrix), ]
# # 确保字符型列是因子（因子便于颜色映射）
# rightdf$Gene_ID <- factor(rightdf$Gene_ID)
# rightdf$SNP_Annotation <- factor(rightdf$SNP_Annotation)
# rightdf$COG_category <- factor(rightdf$COG_category)
# # 定义颜色映射规则（基于因子级别）
# gene_id_colors <- setNames(
#   pal_lancet("lanonc")(length(levels(rightdf$Gene_ID))), # 提取对应数量的颜色
#   levels(rightdf$Gene_ID)
# )
# snp_annotation_colors <- setNames(
#   pal_jama("default")(length(levels(rightdf$SNP_Annotation))), # 提取对应数量的颜色
#   levels(rightdf$SNP_Annotation)
# )
# cog_category_colors <- setNames(
#   pal_npg("nrc")(length(levels(rightdf$COG_category))), # 提取对应数量的颜色
#   levels(rightdf$COG_category)
# )
# # 创建右侧注释条
# right_annotation <- rowAnnotation(
#   df = rightdf, # 注释数据框
#   col = list(
#     Gene_ID = gene_id_colors,          # Gene_ID 的颜色映射
#     SNP_Annotation = snp_annotation_colors,    # SNP_Annotation 的颜色映射
#     COG_category = cog_category_colors  # COG_category 的颜色映射
#   ),
#   show_annotation_name = FALSE
# )
# ht2 <- Heatmap(
#   block_2_matrix,
#   name = "LD R^2", 
#   col = block_colormap, 
#   cluster_rows = FALSE, 
#   cluster_columns = FALSE, 
#   show_row_names = FALSE, 
#   show_column_names = FALSE, 
#   column_dend_height = unit(20, "mm"),
#   show_row_dend = FALSE,
#   left_annotation = left_annotation,
#   right_annotation = right_annotation
# )
# draw(
#   ht2,
#   heatmap_legend_side = "bottom",          
#   annotation_legend_side = "bottom",
#   merge_legend = TRUE
# )
# grid.text(
#   "Block 2",                     # 替换为您的标题
#   x = unit(0.3, "npc"),                   # 标题水平居中
#   y = unit(0.98, "npc"),                  # 标题在顶部，稍微留出间距
#   gp = gpar(fontsize = 16, fontface = "bold") # 字体大小和加粗
# )
# 
# ##
# 
# block_2_0_1matrix<-as.matrix(peakdf[rownames(block_2_matrix),])
# #row_clustering <- hclust(dist(block_2_0_1matrix))
# col_clustering <- hclust(dist(t(block_2_0_1matrix))) 
# block_2_0_1matrix <- apply(block_2_0_1matrix, c(1, 2), function(x) factor(x, levels = c(0, 1)))
# 
# strain_colors <- c("#FF0000", "#800080", "#8A2BE2", "#7FBF7F", "#00FF00")
# variant_colors <- c("0" = "#4575B4", "1" = "#D73027")
# 
# # 优化颜色映射
# color_mapping <- colorRamp2(
#   seq(min(annostrain$strainloading), 
#       max(annostrain$strainloading), 
#       length.out = length(strain_colors)), 
#   strain_colors
# )
# 
# strains_annotation <- HeatmapAnnotation(
#   Strain = annostrain$strainloading,  
#   col = list(Strain = color_mapping),
#   show_annotation_name = FALSE
# )
# 
# ht2.1 <- Heatmap(
#   block_2_0_1matrix,
#   name = "Variants", 
#   cluster_rows = FALSE,  # 使用保存的行聚类树
#   cluster_columns = as.dendrogram(col_clustering),  # 使用保存的列聚类树
#   show_row_names = FALSE, 
#   show_column_names = FALSE, 
#   #show_row_dend = FALSE,
#   left_annotation = left_annotation,
#   right_annotation = right_annotation,
#   top_annotation = strains_annotation,
#   col = variant_colors  # 应用颜色映射
# )
# 
# draw(
#   ht2.1,
#   heatmap_legend_side = "bottom",          
#   annotation_legend_side = "bottom",
#   merge_legend = TRUE
# )
# 
# ###
# #3
# block_3_idx <- which(block_assignment == 3)
# block_3_matrix <- Vr2_matrix[block_3_idx, block_3_idx]
# block_3_matrix <- block_3_matrix[hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="3","label"],
#                                  hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="3","label"]]
# block_3_annotation <- annoall[rownames(block_3_matrix), ] # 根据行索引提取注释数据
# # 定义左侧注释条（rowAnnotation）
# left_annotation <- rowAnnotation(
#   df = block_3_annotation, # 注释数据框
#   col = list(
#     type = type_colors,       # 定义 type 注释的颜色映射
#     PC1 = pc1_colors,         # 定义 PC1 注释的颜色映射
#     position = position_colors # 定义 position 注释的颜色映射
#   ),
#   show_annotation_name = FALSE,
#   annotation_legend_param = list(
#     type = list(title = "Type"),     # 修改 type 图例的标题
#     PC1 = list(title = "PC1^2"), # 修改 PC1 图例的标题
#     position = list(title = "Position") # 修改 position 图例的标题
#   )
# )
# rightdf<-hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block==3,
#                                    c("Gene_ID","SNP_Annotation","COG_category")]
# rightdf <- rightdf[rownames(block_3_matrix), ]
# # 确保字符型列是因子（因子便于颜色映射）
# rightdf$Gene_ID <- factor(rightdf$Gene_ID)
# rightdf$SNP_Annotation <- factor(rightdf$SNP_Annotation)
# rightdf$COG_category <- factor(rightdf$COG_category)
# # 定义颜色映射规则（基于因子级别）
# gene_id_colors <- setNames(
#   pal_lancet("lanonc")(length(levels(rightdf$Gene_ID))), # 提取对应数量的颜色
#   levels(rightdf$Gene_ID)
# )
# snp_annotation_colors <- setNames(
#   pal_jama("default")(length(levels(rightdf$SNP_Annotation))), # 提取对应数量的颜色
#   levels(rightdf$SNP_Annotation)
# )
# cog_category_colors <- setNames(
#   pal_npg("nrc")(length(levels(rightdf$COG_category))), # 提取对应数量的颜色
#   levels(rightdf$COG_category)
# )
# # 创建右侧注释条
# right_annotation <- rowAnnotation(
#   df = rightdf, # 注释数据框
#   col = list(
#     Gene_ID = gene_id_colors,          # Gene_ID 的颜色映射
#     SNP_Annotation = snp_annotation_colors,    # SNP_Annotation 的颜色映射
#     COG_category = cog_category_colors  # COG_category 的颜色映射
#   ),
#   show_annotation_name = FALSE
# )
# ht3 <- Heatmap(
#   block_3_matrix,
#   name = "LD R^2", 
#   col = block_colormap, 
#   cluster_rows = FALSE, 
#   cluster_columns = FALSE, 
#   show_row_names = FALSE, 
#   show_column_names = FALSE, 
#   column_dend_height = unit(20, "mm"),
#   show_row_dend = FALSE,
#   left_annotation = left_annotation,
#   right_annotation = right_annotation
# )
# draw(
#   ht3,
#   heatmap_legend_side = "bottom",          
#   annotation_legend_side = "bottom",
#   merge_legend = TRUE
# )
# grid.text(
#   "Block 3",                     # 替换为您的标题
#   x = unit(0.5, "npc"),                   # 标题水平居中
#   y = unit(0.1, "npc"),                  # 标题在顶部，稍微留出间距
#   gp = gpar(fontsize = 16, fontface = "bold") # 字体大小和加粗
# )
# 
# ##
# 
# block_3_0_1matrix<-as.matrix(peakdf[rownames(block_3_matrix),])
# #row_clustering <- hclust(dist(block_3_0_1matrix))
# col_clustering <- hclust(dist(t(block_3_0_1matrix))) 
# block_3_0_1matrix <- apply(block_3_0_1matrix, c(1, 2), function(x) factor(x, levels = c(0, 1)))
# 
# strain_colors <- c("#FF0000", "#800080", "#8A2BE2", "#7FBF7F", "#00FF00")
# variant_colors <- c("0" = "#4575B4", "1" = "#D73027")
# 
# # 优化颜色映射
# color_mapping <- colorRamp2(
#   seq(min(annostrain$strainloading), 
#       max(annostrain$strainloading), 
#       length.out = length(strain_colors)), 
#   strain_colors
# )
# 
# strains_annotation <- HeatmapAnnotation(
#   Strain = annostrain$strainloading,  
#   col = list(Strain = color_mapping),
#   show_annotation_name = FALSE
# )
# 
# ht3.1 <- Heatmap(
#   block_3_0_1matrix,
#   name = "Variants", 
#   cluster_rows = FALSE,  # 使用保存的行聚类树
#   cluster_columns = as.dendrogram(col_clustering),  # 使用保存的列聚类树
#   show_row_names = FALSE, 
#   show_column_names = FALSE, 
#   show_row_dend = FALSE,
#   left_annotation = left_annotation,
#   right_annotation = right_annotation,
#   top_annotation = strains_annotation,
#   col = variant_colors  # 应用颜色映射
# )
# 
# draw(
#   ht3.1,
#   heatmap_legend_side = "bottom",          
#   annotation_legend_side = "bottom",
#   merge_legend = TRUE
# )
# 
# ###
# #4
# block_4_idx <- which(block_assignment == 4)
# block_4_matrix <- Vr2_matrix[block_4_idx, block_4_idx]
# block_4_matrix <- block_4_matrix[hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="4","label"],
#                                  hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="4","label"]]
# block_4_annotation <- annoall[rownames(block_4_matrix), ] # 根据行索引提取注释数据
# # 定义左侧注释条（rowAnnotation）
# left_annotation <- rowAnnotation(
#   df = block_4_annotation, # 注释数据框
#   col = list(
#     type = type_colors,       # 定义 type 注释的颜色映射
#     PC1 = pc1_colors,         # 定义 PC1 注释的颜色映射
#     position = position_colors # 定义 position 注释的颜色映射
#   ),
#   show_annotation_name = FALSE,
#   annotation_legend_param = list(
#     type = list(title = "Type"),     # 修改 type 图例的标题
#     PC1 = list(title = "PC1^2"), # 修改 PC1 图例的标题
#     position = list(title = "Position") # 修改 position 图例的标题
#   )
# )
# rightdf<-hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block==4,
#                                    c("Gene_ID","SNP_Annotation","COG_category")]
# rightdf <- rightdf[rownames(block_4_matrix), ]
# # 确保字符型列是因子（因子便于颜色映射）
# rightdf$Gene_ID <- factor(rightdf$Gene_ID)
# rightdf$SNP_Annotation <- factor(rightdf$SNP_Annotation)
# rightdf$COG_category <- factor(rightdf$COG_category)
# # 定义颜色映射规则（基于因子级别）
# gene_id_colors <- setNames(
#   pal_lancet("lanonc")(length(levels(rightdf$Gene_ID))), # 提取对应数量的颜色
#   levels(rightdf$Gene_ID)
# )
# snp_annotation_colors <- setNames(
#   pal_jama("default")(length(levels(rightdf$SNP_Annotation))), # 提取对应数量的颜色
#   levels(rightdf$SNP_Annotation)
# )
# cog_category_colors <- setNames(
#   pal_npg("nrc")(length(levels(rightdf$COG_category))), # 提取对应数量的颜色
#   levels(rightdf$COG_category)
# )
# # 创建右侧注释条
# right_annotation <- rowAnnotation(
#   df = rightdf, # 注释数据框
#   col = list(
#     Gene_ID = gene_id_colors,          # Gene_ID 的颜色映射
#     SNP_Annotation = snp_annotation_colors,    # SNP_Annotation 的颜色映射
#     COG_category = cog_category_colors  # COG_category 的颜色映射
#   ),
#   show_annotation_name = FALSE
# )
# ht4 <- Heatmap(
#   block_4_matrix,
#   name = "LD R^2", 
#   col = block_colormap, 
#   cluster_rows = FALSE, 
#   cluster_columns = FALSE, 
#   show_row_names = FALSE, 
#   show_column_names = FALSE, 
#   column_dend_height = unit(20, "mm"),
#   show_row_dend = FALSE,
#   left_annotation = left_annotation,
#   right_annotation = right_annotation
# )
# draw(
#   ht4,
#   heatmap_legend_side = "bottom",          
#   annotation_legend_side = "bottom",
#   merge_legend = TRUE
# )
# grid.text(
#   "Block 4",                     # 替换为您的标题
#   x = unit(0.5, "npc"),                   # 标题水平居中
#   y = unit(0.1, "npc"),                  # 标题在顶部，稍微留出间距
#   gp = gpar(fontsize = 16, fontface = "bold") # 字体大小和加粗
# )
# 
# ##
# 
# block_4_0_1matrix<-as.matrix(peakdf[rownames(block_4_matrix),])
# #row_clustering <- hclust(dist(block_4_0_1matrix))
# col_clustering <- hclust(dist(t(block_4_0_1matrix))) 
# block_4_0_1matrix <- apply(block_4_0_1matrix, c(1, 2), function(x) factor(x, levels = c(0, 1)))
# 
# strain_colors <- c("#FF0000", "#800080", "#8A2BE2", "#7FBF7F", "#00FF00")
# variant_colors <- c("0" = "#4575B4", "1" = "#D73027")
# 
# # 优化颜色映射
# color_mapping <- colorRamp2(
#   seq(min(annostrain$strainloading), 
#       max(annostrain$strainloading), 
#       length.out = length(strain_colors)), 
#   strain_colors
# )
# 
# strains_annotation <- HeatmapAnnotation(
#   Strain = annostrain$strainloading,  
#   col = list(Strain = color_mapping),
#   show_annotation_name = FALSE
# )
# 
# ht4.1 <- Heatmap(
#   block_4_0_1matrix,
#   name = "Variants", 
#   cluster_rows = FALSE,  # 使用保存的行聚类树
#   cluster_columns = as.dendrogram(col_clustering),  # 使用保存的列聚类树
#   show_row_names = FALSE, 
#   show_column_names = FALSE, 
#   #show_row_dend = FALSE,
#   left_annotation = left_annotation,
#   right_annotation = right_annotation,
#   top_annotation = strains_annotation,
#   col = variant_colors  # 应用颜色映射
# )
# 
# draw(
#   ht4.1,
#   heatmap_legend_side = "bottom",          
#   annotation_legend_side = "bottom",
#   merge_legend = TRUE
# )

###
#5
block_5_idx <- which(block_assignment == "Other 40 blocks")
block_5_matrix <- Vr2_matrix[block_5_idx, block_5_idx]
block_5_matrix <- block_5_matrix[hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="Other 40 blocks","label"],
                                 hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="Other 40 blocks","label"]]
block_5_annotation <- annoall[rownames(block_5_matrix), ] # 根据行索引提取注释数据
# 定义左侧注释条（rowAnnotation）
left_annotation <- rowAnnotation(
  df = block_5_annotation, # 注释数据框
  col = list(
    type = type_colors,       # 定义 type 注释的颜色映射
    PC1 = pc1_colors,         # 定义 PC1 注释的颜色映射
    position = position_colors # 定义 position 注释的颜色映射
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    type = list(title = "Type",
                title_gp = gpar(fontsize = 13),  
                labels_gp = gpar(fontsize = 12)),     
    PC1 = list(title = expression(PC1^2),
               at=seq(300,1000,by=300),
               title_gp = gpar(fontsize = 13),  
               labels_gp = gpar(fontsize = 12)), 
    position = list(
      title = "Position (Mb)",  # 修改 position 图例的标题
      at = seq(position_min, position_max, by = 1000000),  # 设置 position 的显示范围和间隔
      labels = c("0", "1", "2", "3", "4", "5"),  # 设置图例标签
      title_gp = gpar(fontsize = 13),  
      labels_gp = gpar(fontsize = 12)) 
  )
)
rightdf<-hitsinfokeep_sortedblocks[hitsinfokeep_sortedblocks$Block=="Other 40 blocks",
                                   c("Gene_ID","SNP_Annotation","COG_category")]
rightdf <- rightdf[rownames(block_5_matrix), ]

rightdf$Blocks<-factor(c(rep(1,9),rep(2,5),3,4,4,rep(5,4),rep(6,5),rep(7,7),c(8:13),14,14,15,rep(16,6),
                 17,17,18,18,18,19,19,19,rep(20,7),21,21,22,23,24,25,rep(26,4),27,27,c(28:40)))
r<-rownames(rightdf)
b3<-hitsanno_sorted[hitsanno_sorted$Block == "other 40 blocks", ]
library(dplyr)
b3 <- b3 %>%
  mutate(
    `old locus tag(NC_012731.1)` = case_when(
      Gene_ID == "KP1_RS30365" ~ Gene_ID,
      Gene_ID == "KP1_RS30360" ~ Gene_ID,
      TRUE ~ `old locus tag(NC_012731.1)`  # 保持原值
    )
  )
rightdf$Gene_ID<-b3$`old locus tag(NC_012731.1)`
rownames(rightdf)==r
rightdf <- rightdf[rownames(block_5_matrix), ]
# 确保字符型列是因子（因子便于颜色映射）
rightdf$Gene_ID <- factor(rightdf$Gene_ID)
rightdf$SNP_Annotation <- factor(rightdf$SNP_Annotation)
rightdf$COG_category <- factor(rightdf$COG_category)
# 定义颜色映射规则（基于因子级别）
gene_id_colors <- setNames(
  colorRampPalette(pal_lancet("lanonc")(9))(length(levels(rightdf$Gene_ID))), # 提取对应数量的颜色
  levels(rightdf$Gene_ID)
)
snp_annotation_colors <- setNames(
  pal_jama("default")(length(levels(rightdf$SNP_Annotation))), # 提取对应数量的颜色
  levels(rightdf$SNP_Annotation)
)
cog_category_colors <- setNames(
  colorRampPalette(pal_npg("nrc")(10))(length(levels(rightdf$COG_category))), # 提取对应数量的颜色
  levels(rightdf$COG_category)
)
# 创建右侧注释条
right_annotation <- rowAnnotation(
  df = rightdf[,1:3], # 注释数据框
  col = list(
    Gene_ID = gene_id_colors,          # Gene_ID 的颜色映射
    SNP_Annotation = snp_annotation_colors,    # SNP_Annotation 的颜色映射
    COG_category = cog_category_colors  # COG_category 的颜色映射
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Gene_ID = list(title = "Gene",
                   title_gp = gpar(fontsize = 13),  
                   labels_gp = gpar(fontsize = 12),
                   ncol = 5),     
    SNP_Annotation = list(title = "SNP",
                          title_gp = gpar(fontsize = 13),  
                          labels_gp = gpar(fontsize = 12)), 
    COG_category = list(
      title = "COG",  
      title_gp = gpar(fontsize = 13),  
      labels_gp = gpar(fontsize = 12),
      ncol = 4) 
  )
)
library(RColorBrewer)  
palette1 <- brewer.pal(8, "Set3")  # 调色板1
palette2 <- brewer.pal(8, "Dark2")  # 调色板2
palette3 <- brewer.pal(8, "Paired")  # 调色板3
palette4 <- brewer.pal(8, "Accent")  # 调色板4
palette5 <- brewer.pal(8, "Pastel1")  # 调色板5
# 合并调色板，生成40种颜色
all_colors <- c(palette1, palette2, palette3, palette4, palette5)
# 为每个block分配一种颜色
block_colors <- setNames(
  all_colors[1:length(unique(rightdf$Blocks))],  # 只取前40种颜色
  as.character(unique(rightdf$Blocks))
)
ha_right <- HeatmapAnnotation(
  Blocks = rightdf$Blocks,
  col = list(Blocks = block_colors),
  show_annotation_name = TRUE,
  show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 13)
)

ht5 <- Heatmap(
  block_5_matrix,
  name = "LD R^2", 
  col = block_colormap, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  column_dend_height = unit(20, "mm"),
  show_row_dend = FALSE,
  left_annotation = left_annotation,
  top_annotation = ha_right,
  #column_split = 40,
  #row_split = 40
  right_annotation = right_annotation,
  heatmap_legend_param = list(
    title = expression("LD " * R^2),  # 设置图例标题
    title_gp = gpar(fontsize = 13),  
    labels_gp = gpar(fontsize = 12) 
  )
)###40
draw(
  ht5,
  heatmap_legend_side = "bottom",          
  annotation_legend_side = "bottom",
  merge_legend = TRUE,
  padding = unit(c(0.2, 0, 1, 0), "cm") 
)
grid.text(
  "Other 40 blocks",                     # 替换为您的标题
  x = unit(0.5, "npc"),                   # 标题水平居中
  y = unit(0.97, "npc"),                  # 标题在顶部，稍微留出间距
  gp = gpar(fontsize = 16, fontface = "bold") # 字体大小和加粗
)

##
pcvalues<-pca[order(pca$PC1),"ID"]
block_5_0_1matrix<-as.matrix(peakdf[rownames(block_5_matrix),pcvalues])
#row_clustering <- hclust(dist(block_5_0_1matrix))
#col_clustering <- hclust(dist(t(block_5_0_1matrix))) 
block_5_0_1matrix <- apply(block_5_0_1matrix, c(1, 2), function(x) factor(x, levels = c(0, 1)))

#strain_colors <- c("#FF0000", "#800080", "#8A2BE2", "#7FBF7F", "#00FF00")
variant_colors <- c("0" = "#4575B4", "1" = "#D73027")

# 优化颜色映射
color_mapping <- circlize::colorRamp2(
  seq(min(annostrain$strainloading), max(annostrain$strainloading), length.out = 369),  # 创建渐变范围
  viridis(369, option = "H")  # 使用 viridis 的 "H" 配色
)

strains_annotation <- HeatmapAnnotation(
  Strain = annostrain$strainloading,  
  col = list(Strain = color_mapping),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Strain = list(title = "Strain",
                  at=seq(-0.3,0.4,by=0.2),
                  title_gp = gpar(fontsize = 13),  
                  labels_gp = gpar(fontsize = 12))
  )
)

ht5.1 <- Heatmap(
  block_5_0_1matrix,
  name = "Variants", 
  cluster_rows = FALSE,  
  cluster_columns = FALSE,  
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  show_row_dend = FALSE,
  left_annotation = left_annotation,
  right_annotation = right_annotation,
  top_annotation = strains_annotation,
  col = variant_colors,  # 应用颜色映射
  heatmap_legend_param = list(
    title = "Variant",  # 设置图例标题
    title_gp = gpar(fontsize = 13),  
    labels_gp = gpar(fontsize = 12) 
  )
)

draw(
  ht5.1,
  heatmap_legend_side = "bottom",          
  annotation_legend_side = "bottom",
  merge_legend = TRUE,
  padding = unit(c(0.2, 0, 1, 0), "cm") 
)
grid.text(
  "Other 40 blocks",                     # 替换为您的标题
  x = unit(0.5, "npc"),                   # 标题水平居中
  y = unit(0.97, "npc"),                  # 标题在顶部，稍微留出间距
  gp = gpar(fontsize = 16, fontface = "bold") # 字体大小和加粗
)

###
z1<-ggplot(squdf,aes(x=POS,y=PC1)) +
  geom_point(aes(color = PC1 > 300), alpha = 0.5, size = 3) +  # 根据条件设置颜色
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept = 300, color = "blue", linetype = "dashed", linewidth = 1) +  # 添加水平线，y=250
  labs(x = "Position (bp)", y = "Peak 1") +
  theme_classic()+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")+
  xlim(343000,353000)+
  geom_text_repel(data = x,aes(label=old_locus_tag),size=5)+
  scale_y_continuous(breaks = seq(0, 1000, by = 300))

z2<-ggplot(squdf,aes(x=POS,y=PC1)) +
  geom_point(aes(color = PC1 > 300), alpha = 0.5, size = 3) +  # 根据条件设置颜色
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept = 300, color = "blue", linetype = "dashed", linewidth = 1) +  # 添加水平线，y=250
  labs(x = "Position (bp)", y = "Peak 2") +
  theme_classic()+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")+
  geom_text_repel(data = x,aes(label=old_locus_tag),size=5)+
  xlim(1568000,1578000)+
  scale_y_continuous(breaks = seq(0, 1000, by = 300))

z3<-ggplot(squdf,aes(x=POS,y=PC1)) +
  geom_point(aes(color = PC1 > 300), alpha = 0.5, size = 3) +  # 根据条件设置颜色
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept = 300, color = "blue", linetype = "dashed", linewidth = 1) +  # 添加水平线，y=250
  labs(x = "Position (bp)", y = "Peak 3") +
  theme_classic()+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")+
  xlim(1928000,1938000)+
  geom_text_repel(data = x,aes(label=old_locus_tag),size=5)+
  scale_y_continuous(breaks = seq(0, 1000, by = 300))

z4<-ggplot(squdf,aes(x=POS,y=PC1)) +
  geom_point(aes(color = PC1 > 300), alpha = 0.5, size = 3) +  # 根据条件设置颜色
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept = 300, color = "blue", linetype = "dashed", linewidth = 1) +  # 添加水平线，y=250
  labs(x = "Position (bp)", y = "Peak 4") +
  theme_classic()+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")+
  xlim(5125000,5135000)+
  geom_text_repel(data = x,aes(label=old_locus_tag),size=5)+
  scale_y_continuous(breaks = seq(0, 1000, by = 300))

library(cowplot)
plot_grid(z1,z2,z3,z4,
          align = "hv",ncol = 1)

###heatmap
infodf<-squdf[,c("POS","PC1")]
infodf<-infodf[infodf$PC1>300,]
hm<-geno[rownames(infodf),]
library(pheatmap)
pheatmap(as.matrix(hm),
         show_rownames = FALSE,show_colnames = FALSE)
###50k
squtable<-squdf[c("PC1","POS")]
squtable$maxPC1_50k<-NA
for (i in 1:nrow(squtable)) {
  temp<-squtable[(squtable$POS<(squtable[i,2]+50000))&(squtable$POS>=squtable[i,2]),]
  squtable[i,3]<-max(temp$PC1)
}

hist(squtable$maxPC1_50k,col='blue',border='yellow',
     main='Highest loading in the vicinity, window size = 50kb',
     xlab='Squared PC1 loading values',breaks = 50,
     cex.main=1.5,     # 控制主标题的字体大小
     cex.lab=1.2,      # 控制x轴和y轴标签的字体大小
     cex.axis=1.1)
abline(v=c(160, 240), col="red", lwd=2) 

ggplot(squtable, aes(x = maxPC1_50k)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "blue", alpha = 0.5, color = "yellow") +  # 直方图
  geom_density(color = "blue", size = 1) +  # 密度曲线
  geom_vline(xintercept = c(160, 240), color = "red", linetype = "dashed", size = 1) +  # 添加参考线
  labs(
    title = "Highest loading in the vicinity, window size = 50kb",
    x = "Squared PC1 loading values",
    y = "Density"
  ) +
  theme_minimal() +  # 简洁主题
  theme(#axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 12),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 12),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 14),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 14), # 设置y轴标题的字体大小
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none")+
  # 添加标注
  annotate("text", x = 20, y = 0.006, label = "Low", color = "black", size = 8, hjust = 0) +
  annotate("text", x = 170, y = 0.007, label = "Middle", color = "black", size = 4, hjust = 0) +
  annotate("text", x = 550, y = 0.004, label = "High", color = "black", size = 8, hjust = 0)

# arrows(70, 3000, 30, 6000, col="blue", lwd=2, length=0.1)  # 调整箭头位置
# text(70, 2000, "Low loading SNPs", col="blue", cex=1.5, font=2, pos=2)  # 文本位置调整
# arrows(350, 4000, 200, 8000, col="blue", lwd=2, length=0.1)  # 调整箭头位置
# text(350, 5000, "High loading SNPs", col="blue", cex=1.5, font=2, pos=4)  # 文本位置调整

h<-squtable[squtable$maxPC1_50k>240,]
h<-h[c("POS")]
#write.table(h,file = "half_matching/halfloci/highloadings_pos.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)
l<-squtable[squtable$maxPC1_50k<160,]
l<-l[c("POS")]
#write.table(l,file = "half_matching/halfloci/lowloadings_pos.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)

###anno snp
anno1<-read.csv("NR369maf002anno.txt",sep = "|",header = TRUE,fill = TRUE)
anno2<-read.csv("MM_yezscq0m.emapper.annotations.tsv",sep = "\t")
library(dplyr)
snp_info <-left_join(anno1, anno2)
squdf_info<-left_join(squdf,snp_info)
x<-squdf_info[squdf_info$PC1>300,]
x<-x[((x$POS>300000)&(x$POS<400000))|((x$POS>1520000)&(x$POS<1620000))|((x$POS>1880000)&(x$POS<1980000))|((x$POS>5080000)&(x$POS<5180000)),]
# x<-x[x$Annotation!="synonymous_variant",]
# x<-x[x$Gene_ID!="null",]
result <- x %>%
  group_by(Gene_ID) %>%
  summarize(PC1 = max(PC1),
            pos_count = n())

numberofnonsys<-x[,c("Gene_ID","Annotation")]
numberofnonsys<-numberofnonsys[numberofnonsys$Annotation!="synonymous_variant",]
table(numberofnonsys$Gene_ID)

x<-left_join(result,x)
x<-x[order(x$POS),]
x<-x[!duplicated(x[,c("PC1","Gene_ID")]),]
x$geneName <- gsub("gene-", "", x$Gene_ID)
x <- x[-c(3), ]
library(dplyr)
x <- x %>%
  mutate(
    old_locus_tag = case_when(
      geneName == "KP1_RS30365" ~ geneName,
      geneName == "KP1_RS30360" ~ geneName,
      TRUE ~ old_locus_tag  # 保持原值
    )
  )

t<-x[,c("Gene_ID","POS","gene","product","COG_category")]

library(ggrepel)
ggplot(load_snp_df,aes(x=POS,y=PC1))+
  geom_point()+
  theme_classic()+
  xlab("Position (Mb)")+
  ylab("PC1 Loading")+
  scale_x_continuous(labels = c("0","1","2","3","4","5"))+
  geom_text_repel(data = x,aes(label=old_locus_tag),max.overlaps=20)+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16))
###
anno1<-read.csv("NR369maf002anno.txt",sep = "|",header = TRUE,fill = TRUE)
anno2<-read.csv("MM_yezscq0m.emapper.annotations.tsv",sep = "\t")
library(dplyr)
snp_info <-left_join(anno1, anno2)
squdf_info<-left_join(load_snp_df,snp_info)
coresnp<-squdf_info[,c("POS","PC1","Gene_ID","old_locus_tag","COG_category")]
coresnp <- coresnp %>%
  mutate(
    old_locus_tag = ifelse(
      is.na(old_locus_tag) | old_locus_tag == "",
      gsub("^gene-", "", Gene_ID),
      old_locus_tag
    )
  )
coresnp_clean <- coresnp %>%
  # 1. 去掉 NA 和 "-" 的 COG_category
  filter(!is.na(COG_category), COG_category != "-") %>%
  
  # 2. 拆分成单个 COG 字母
  mutate(COG_category = strsplit(COG_category, "")) %>%
  unnest(COG_category) %>%
  
  # 3. 转换 PC1 为绝对值
  mutate(PC1 = abs(PC1))

mean_df <- coresnp_clean %>%
  group_by(old_locus_tag) %>%
  summarise(avg_PC1 = mean(PC1, na.rm = TRUE), .groups = "drop")

# 2. 每个 old_locus_tag 保留一个代表行（含 Gene_ID、COG_category）
gene_info <- coresnp_clean %>%
  group_by(old_locus_tag) %>%
  select(old_locus_tag, Gene_ID, COG_category)%>%
  distinct()

# 3. 合并
meangene <- right_join(mean_df, gene_info, by = "old_locus_tag")
meangene <- meangene %>%
  left_join(cog_function, by = "COG_category")

summary_df <- meangene %>%
  group_by(COG_category) %>%
  summarise(
    mean_PC1 = mean(avg_PC1, na.rm = TRUE),
    sd_PC1   = sd(avg_PC1, na.rm = TRUE),
    se_PC1sq   = sd(avg_PC1, na.rm = TRUE) / sqrt(n()), 
    n        = n(),
    Category = first(Category),  # 或用 unique(Category)[1]
    Function = first(Function),  # 你要是想按 Function 展示位置
    .groups = "drop"
  )

summary_df <- summary_df %>%
  arrange(Category, desc(mean_PC1)) %>%
  mutate(Function = factor(Function, levels = unique(Function)))
meangene <- meangene %>%
  mutate(Function = factor(Function, levels = levels(summary_df$Function)))

psnp<-ggplot() +
  geom_jitter(
    data = meangene,
    aes(x = Function, y = avg_PC1, color = Category),
    width = 0.2, alpha = 0.5, size = 2, show.legend = FALSE
  ) +
  geom_errorbar(
    data = summary_df,
    aes(
      x = Function,
      ymin = mean_PC1 - se_PC1sq,
      ymax = mean_PC1 + se_PC1sq
    ),
    width = 0.3, linewidth = 0.6
  ) +
  # 平均值横线
  geom_errorbar(
    data = summary_df,
    aes(
      x = Function,
      ymin = mean_PC1,
      ymax = mean_PC1
    ),
    width = 0.5, linewidth = 1
  )+
  coord_flip() +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_color_npg() +  
  labs(
    x = NULL,
    y = "Absolute Value of PC1 Loadings",
    title = "Core Genes"
  ) +
  theme_classic() +
  theme(
    strip.text.y = element_text(size = 12, angle = 0, hjust = 0),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5)
  )+
  scale_x_discrete(limits = rev)+
  geom_hline(yintercept = mean(meangene$avg_PC1),linewidth = 0.5,linetype = "solid",color="#EE0000FF")+
  scale_y_continuous(limits = c(0, max(geneloading_anno$PC1)),
                     breaks = seq(0, max(geneloading_anno$PC1),10))

plot_grid(psnp,pgene,align = "hv",ncol = 1)
###
my_mean_se <- function(x) {
  m <- mean(x)
  s <- sd(x)
  n <- length(x)
  se <- s / sqrt(n)
  return(c(y = m, ymin = m - se, ymax = m + se))
}
library(ggpubr)
ggplot(meangene, aes(x = Category, y = avg_PC1, fill = Category)) +
  stat_summary(fun = mean, geom = "bar", width = 0.5, alpha = 0.6, color = "black") +
  stat_summary(fun.data = my_mean_se, geom = "errorbar", width = 0.2, color = "black") +
  geom_jitter(width = 0.1, size = 2, alpha = 0.1, color = "black") +
  stat_compare_means(
    method = "t.test",
    comparisons = list(
      c("CELLULAR PROCESSES AND SIGNALING", "INFORMATION STORAGE AND PROCESSING"),
      c("CELLULAR PROCESSES AND SIGNALING", "METABOLISM"),
      c("INFORMATION STORAGE AND PROCESSING", "METABOLISM"),
      c("INFORMATION STORAGE AND PROCESSING", "POORLY CHARACTERIZED"),
      c("CELLULAR PROCESSES AND SIGNALING", "POORLY CHARACTERIZED"),
      c("METABOLISM", "POORLY CHARACTERIZED")
    ),
    label = "p.format"
  ) +
  theme_light() +
  labs(y = "Projected PC1", x = NULL) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = "none"
  )+scale_fill_npg()
####
peakinfo<-squdf_info[peak_positions,]
peakinfo<-peakinfo[!((peakinfo$POS>5124000)&(peakinfo$POS<5130000)|
                       (peakinfo$POS>1929000)&(peakinfo$POS<1937000)|
                       (peakinfo$POS>1571000)&(peakinfo$POS<1578000)|
                       (peakinfo$POS>346000)&(peakinfo$POS<353000)),]
peakinfo<-peakinfo[peakinfo$Gene_ID!="null",]
t2<-peakinfo[,c("Gene_ID","POS","gene","product","COG_category","Annotation")]
t2$Gene_ID <- gsub("gene-", "", t2$Gene_ID)
write.csv(t2,"figs/suptable.csv",
            row.names = FALSE,quote = FALSE)
# ##
# loaddf_info<-left_join(load_snp_df,snp_info)
# x<-loaddf_info[(loaddf_info$PC1>20)|(loaddf_info$PC1< -15),]
# x<-x[x$Gene_ID!="null",]
# result <- x %>%
#   group_by(Gene_ID) %>%
#   summarize(PC1 = max(abs(PC1)))
###mean snp
snpdf_info<-left_join(load_snp_df,snp_info)
dfinfo<-na.omit(snpdf_info[,c("Gene_ID","POS","COG_category","PC1")])

snpdf<-dfinfo[dfinfo$COG_category!="-",]
snpdf <- snpdf %>%
  mutate(COG_category = str_split(COG_category,"")) %>%
  unnest(COG_category)
snpdf$PC1<-abs(snpdf$PC1)

avergeneload<-snpdf %>%
  group_by(Gene_ID) %>%
  summarise(snpcount = n(),
            genemean = mean(PC1),
            COG_category = first(COG_category),
            pos=mean(POS))

snpsum <- avergeneload %>%
  group_by(COG_category) %>%
  summarise(mean = mean(genemean),
            count = n(),
            se = sd(genemean)/sqrt(n()))

snpsum<-snpsum[order(snpsum$mean),]
snpsum$COG_category<-factor(snpsum$COG_category,levels = unique(snpsum$COG_category))

library(ggplot2)
ggplot(snpsum,aes(x=COG_category))+
  geom_bar(aes(y=mean),stat = "identity",fill="lightblue",alpha=0.7)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.2)+
  geom_text(aes(y=mean,label = count),vjust= 3,color="black")+
  labs(title = "COG category Mean, Count, and Standard Error",y = "Mean Absolute Loading")+
  theme_minimal()+
  geom_hline(yintercept = mean(avergeneload$genemean),linetype="dashed",color="red",linewidth=1)

###clusterprofiler
library(clusterProfiler)
library(tidyverse)
gene_cog <- anno2[,c("protein_id","COG_category")]
gene_cog <- gene_cog[gene_cog$COG_category!="-",]
gene_cog <- gene_cog %>%
  mutate(COG_category = str_split(COG_category,"")) %>%
  unnest(COG_category)

target_genes <- unique(x$protein_id)
target_genes <- setdiff(target_genes, "")

cog2gene_list <- split(gene_cog$protein_id, gene_cog$COG_category)
cog2gene <- stack(cog2gene_list)
colnames(cog2gene) <- c("GENE", "TERM")
cog2gene<-cog2gene[,c("TERM","GENE")]

ego <- enricher(
  target_genes,
  TERM2GENE = cog2gene
)
barplot(ego)
dotplot(ego)

###
gene_go <- anno2[,c("protein_id","GOs")]
gene_go <- gene_go[gene_go$GOs!="-",]
gene_go <- gene_go %>%
  mutate(GOs = str_split(GOs,",")) %>%
  unnest(GOs)
cog2go_list <- split(gene_go$protein_id, gene_go$GOs)
cog2go <- stack(cog2go_list)
colnames(cog2go) <- c("GENE", "TERM")
cog2go<-cog2go[,c("TERM","GENE")]
ego2 <- enricher(
  target_genes,
  TERM2GENE = cog2go
)
barplot(ego2)
###
gene_go <- anno2[,c("protein_id","KEGG_ko")]
gene_go <- gene_go[gene_go$KEGG_ko!="-",]
gene_go <- gene_go %>%
  mutate(KEGG_ko = str_split(KEGG_ko,",")) %>%
  unnest(KEGG_ko)
cog2go_list <- split(gene_go$protein_id, gene_go$KEGG_ko)
cog2go <- stack(cog2go_list)
colnames(cog2go) <- c("GENE", "TERM")
cog2go<-cog2go[,c("TERM","GENE")]
ego3 <- enricher(
  target_genes,
  TERM2GENE = cog2go
)
barplot(ego3)
###
y<-squdf_info[squdf_info$PC2>75,]
y<-y[y$Gene_ID!="null",]
result2 <- y %>%
  group_by(Gene_ID) %>%
  summarize(PC2 = max(PC2))
y<-left_join(result2,y)
y<-y[order(y$POS),]
y<-y[!duplicated(y[,c("PC2","Gene_ID")]),]
y$geneName <- gsub("gene-", "", y$Gene_ID)

ggplot(squdf,aes(x=POS,y=PC2,color=POS))+
  geom_point()+
  theme_classic()+
  xlab("Position (Mb)")+
  ylab("Squared PC2 Loading")+
  scale_x_continuous(labels = c("0","1","2","3","4","5"))+
  scale_color_gradientn(colours = c("firebrick", "springgreen4", 
                                    "blue3", "turquoise3","darkorchid2", "gold2","orange"),
                        guide="none")+
  geom_text_repel(data = y,aes(label=geneName),max.overlaps=20)+
  theme(axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16))
###vistree
library(ggtree)
library(treeio)
allKp<-read.newick("KP_1421_accessions.clean.snp.Fast.tree_nt")
metainfo<-read.csv("allKpmeta.csv")
countryinfo <- strsplit(as.character(metainfo$geo_loc_name), ":")
countryinfo1 <- sapply(countryinfo, `[`, 1)
metainfo$country<-countryinfo1
library(countrycode)
metainfo$Continent<-countrycode(metainfo$country, "country.name", "continent")
metainfo$Continent <- ifelse(is.na(metainfo$Continent), "NA", metainfo$Continent)
metainfo[!is.na(metainfo$country) & metainfo$country == "United Kingdom (England, Wales & N. Ireland)", "Continent"] <- "Europe"

y<-metainfo[,c("ID","Continent","virulence_score","resistance_score")]
y$Continent <- ifelse(is.na(y$Continent), "missing", y$Continent)
#write.csv(y,file = "fs/KP1422info.csv",quote = FALSE,row.names = FALSE)

ContinentInfo<-list(Asia=c(metainfo[metainfo$Continent=="Asia","ID"]),
                    Europe=c(metainfo[metainfo$Continent=="Europe","ID"]),
                    `NA`=c(metainfo[metainfo$Continent=="NA","ID"]),
                    America=c(metainfo[metainfo$Continent=="Americas","ID"]),
                    Africa=c(metainfo[metainfo$Continent=="Africa","ID"]),
                    Oceania=c(metainfo[metainfo$Continent=="Oceania","ID"]))

treeGlobal<-groupOTU(allKp,ContinentInfo,group_name = "Continent")

ggtree(treeGlobal,layout = "ape")+
  geom_tippoint(aes(color = Continent, shape = Continent),size=3)+
  scale_color_npg()+
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(1, "lines"),
        legend.position = "bottom")+
  geom_treescale(offset = 0.001,fontsize=5, linesize=1, x=0, y=0,width=0.01)

###nr369
NR369tree<-read.newick("NR369_accessions.clean.snp.Fast.tree_nt")

pca <- read.table("NR369pca.eigenvec")
rownames(pca)<-pca$V1
pca <- pca[,-c(1,2)]
colnames(pca) <- paste0("PC", 1:ncol(pca))
val<-read.table("NR369pca.eigenval")
colnames(val)<-"val"
rownames(val)<-paste0("PC", 1:nrow(val))
pca<-as.data.frame(t(t(pca)*sqrt(val$val)))

pca$label<-rownames(pca)

tree.a<-full_join(NR369tree,pca,by="label")
library(viridis)
fig1a<-ggtree(tree.a, layout = "ape") +  
  geom_tippoint(aes(color = PC1),size=3) +  
  scale_color_viridis(
    option = "H",  # 使用 Viridis 渐变色
    name = "PC1",  # 图例标题
    breaks = seq(-0.3, 0.4, by = 0.2),  # 设置图例断点
    guide = guide_colorbar(reverse = TRUE)
  )+
  geom_treescale(offset = 0.001,fontsize=5, linesize=1, x=0, y=0,width=0.01)+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(2, "lines"),
        legend.position = "right")
ggsave("/Users/47liu/Documents/paper/peer review/figs/figure1a.pdf",
       plot = fig1a,
       width = 8, height = 6,
       units = "in",
       dpi = 600,
       bg = "transparent",
       useDingbats = FALSE  
)
###
blockstree<-read.newick("block1.tree")
tree.b<-full_join(blockstree,pca,by="label")
library(viridis)
ggtree(tree.b, layout = "ape") +  
  geom_tippoint(aes(color = PC1),size=3) +  
  scale_color_viridis(
    option = "H",  # 使用 Viridis 渐变色
    name = "PC1",  # 图例标题
    breaks = seq(-0.3, 0.4, by = 0.2),  # 设置图例断点
    guide = guide_colorbar(reverse = TRUE)
  )+
  geom_treescale(offset = 0.001,fontsize=5, linesize=1, x=0, y=0,width=0.01)+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(2, "lines"),
        legend.position = "left")
###
groupinfo<-read.csv("fs/groups.csv")
groupinfo<-split(groupinfo$label, groupinfo$group)

treeinfo<-groupOTU(NR369tree,groupinfo,group_name = "Groups")

ggtree(treeinfo,aes(color=Groups),layout = "ape")+
  geom_tippoint(aes(shape = Groups),size=3)+
  scale_color_aaas()+
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(1, "lines"),
        legend.position = "bottom")+
  geom_treescale(offset = 0.001,fontsize=5, linesize=1, x=0, y=0,width=0.01)
###
library(ape)
library(ggtree)
library(dplyr)
# 加载你的树
tree <- tree.a  # 假设你已经加载树文件到 tree.a 中
# 绘制 unrooted 树
p <- ggtree(tree, layout = "unrooted")
# 获取节点数据
node_data <- p$data
# 计算节点相对于树中心的角度
# atan2 计算角度（弧度），转换为度数后以 180° 为分界
node_data <- node_data %>%
  mutate(angle = atan2(y, x) * 180 / pi,  # 计算角度
         angle = ifelse(angle < 0, angle + 360, angle),  # 转换到 0-360 度
         group = ifelse(angle >= 114 & angle < 300, "Left", "Right"))  # 划分左右

# 可视化分组结果
ggtree(tree, layout = "unrooted") +
  geom_point(aes(x = x, y = y, color = group), size = 3, data = node_data) +
  scale_color_manual(values = c("Left" = "blue", "Right" = "red")) +
  theme(legend.position = "bottom")
##
leaf_nodes <- node_data %>% filter(isTip)
# 按组统计菌株名称
grouped_strains <- leaf_nodes %>%
  group_by(group) %>%
  summarise(
    count = n(),                     # 统计菌株数量
    strains = paste(label, collapse = ", ")  # 将菌株名称拼接为字符串
  )
blue_strains <- leaf_nodes %>% filter(group == "Left") %>% pull(label)
red_strains <- leaf_nodes %>% filter(group == "Right") %>% pull(label)
# write.table(data.frame(blue_strains), "./half_matching/halfsample/Half184.list", 
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(data.frame(red_strains), "./half_matching/halfsample/Half185.list", 
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
####LD vs decay
load_snp_df <- as.data.frame(t(load_snp))
load_snp_df$POS<-as.numeric(rownames(load_snp_df))
df<-load_snp_df[,c("POS","PC1")]

my_data <- data.frame(
  SNP_distance = seq(500, 30000, by = 500),
  Pearson_Correlation_Coefficient = numeric(60)
)

for (i in 1:60) {
  print(i)
  p2_loading <- numeric(nrow(df))
  for (j in 1:nrow(df)) {
    p1 <- df$POS[j]
    p2_range <- df$POS >= p1 + 500 * (i-1) & df$POS <= p1 + 500 * i
    p2_loading[j] <- mean(df$PC1[p2_range])
  }
  my_data$Pearson_Correlation_Coefficient[i] <- cor(df$PC1, p2_loading,use = "pairwise.complete.obs")
}
####
window_size <- 100000 
num_chunks <- ceiling(max(df$POS) / window_size)
my_data3 <- data.frame(
  SNP_window = seq(window_size, num_chunks*window_size, by = window_size),
  sumloading = numeric(num_chunks)
)
for (i in 1:num_chunks) {
  print(i)
  start_row <- (i - 1) * window_size + 1  
  end_row <- min(i * window_size, max(df$POS))
  chunk_range <- df$POS >= start_row & df$POS <= end_row
  my_data3$sumloading[i] <- sum(df$PC1[chunk_range])
}
plot(my_data3$SNP_window,my_data3$sumloading,type = "l",
     main = paste0("window size ",window_size/1000," kb"),
     xlab = "SNP window",ylab = "sum(loading)")
###LD
ldws10<-read.table("NR369maf002.matrix.MajorAMinorT.statRR.avg.matrix.step10.All.stat")
ldws500<-read.table("NR369maf002.matrix.MajorAMinorT.statRR.avg.matrix.step500.All.stat")

####strain
calculate_row_frequency <- function(row) {
  num_ones <- sum(row == 1)
  num_zeros <- sum(row == 0)
  return(ifelse(row == 1, num_ones/(num_ones+num_zeros), -(num_zeros/(num_zeros+num_ones))))
}

snpfreq <- as.data.frame(apply(geno, 1, calculate_row_frequency))

snpfreq_normalized<-scale(as.matrix(snpfreq))

snpfreq_normalized <- as.data.frame(snpfreq_normalized)
snpfreq_loading<-t(snpfreq_normalized)*load_snp_df$PC1
snpfreq_loading<-as.data.frame(snpfreq_loading)

r_x<-pca[order(pca$PC1),]
snpfreq_loading<-snpfreq_loading[,rownames(r_x)]
snpfreq_loading$POS<-as.numeric(rownames(snpfreq_loading))

plot(snpfreq_loading$POS,snpfreq_loading$GCF_900084635_1,type = "p")

minstrain<-ggplot(snpfreq_loading,aes(x=POS,y=GCF_900084635_1))+
  geom_point(color="firebrick1")+
  theme_bw()+
  xlab("Position (Mb)")+
  ylab("PC1 loadings")+
  scale_x_continuous(labels = c("0","1","2","3","4","5"))+
  ggtitle("Minimum PC1 strain")+
  ylim(-40,50)

midstrain<-ggplot(snpfreq_loading,aes(x=POS,y=GCF_002187405_1))+
  geom_point(color="yellow1")+
  theme_bw()+
  xlab("Position (Mb)")+
  ylab("PC1 loadings")+
  scale_x_continuous(labels = c("0","1","2","3","4","5"))+
  ggtitle("Middle PC1 strain")+
  ylim(-40,50)

maxstrain<-ggplot(snpfreq_loading,aes(x=POS,y=GCF_016762175_1))+
  geom_point(color="blue1")+
  theme_bw()+
  xlab("Position (Mb)")+
  ylab("PC1 loadings")+
  scale_x_continuous(labels = c("0","1","2","3","4","5"))+
  ggtitle("Maximum PC1 strain")+
  ylim(-40,50)

plot_grid(minstrain,midstrain,maxstrain,align = "hv",ncol = 1)
###decay
my_data2<-data.frame(
  SNP_distance=seq(500, 30000, by = 500),
  maxPC1strain=numeric(60),
  minPC1strain=numeric(60),
  midPC1strain=numeric(60)
)

df2<-snpfreq_loading[,c("POS",
                     "GCF_016762175_1",
                     "GCF_900084635_1",
                     "GCF_002187405_1")]
for (i in 1:60) {
  print(i)
  maxp2_loading<-numeric(nrow(df2))
  minp2_loading<-numeric(nrow(df2))
  midp2_loading<-numeric(nrow(df2))
  for (j in 1:nrow(df2)) {
    p1 <- df2$POS[j]
    p2_range <- df2$POS >= p1 + 500 * (i-1) & df2$POS <= p1 + 500 * i
    maxp2_loading[j] <- mean(df2$GCF_016762175_1[p2_range])
    minp2_loading[j] <- mean(df2$GCF_900084635_1[p2_range])
    midp2_loading[j] <- mean(df2$GCF_002187405_1[p2_range])
  }
  my_data2$maxPC1strain[i] <- cor((df2$GCF_016762175_1), maxp2_loading,use = "pairwise.complete.obs")
  my_data2$minPC1strain[i] <- cor((df2$GCF_900084635_1), minp2_loading,use = "pairwise.complete.obs")
  my_data2$midPC1strain[i] <- cor((df2$GCF_002187405_1), midp2_loading,use = "pairwise.complete.obs")
}
####
window_size <- 10000
num_chunks <- ceiling(max(df2$POS) / window_size)
my_data4 <- data.frame(
  SNP_window = seq(window_size, num_chunks*window_size, by = window_size),
  summaxstrainloading = numeric(num_chunks),
  summinstrainloading = numeric(num_chunks),
  summidstrainloading = numeric(num_chunks)
)
for (i in 1:num_chunks) {
  print(i)
  start_row <- (i - 1) * window_size + 1  
  end_row <- min(i * window_size, max(df$POS))
  chunk_range <- df2$POS >= start_row & df2$POS <= end_row
  my_data4$summaxloading[i] <- sum(df2$GCF_016762175_1[chunk_range])
  my_data4$summinloading[i] <- sum(df2$GCF_900084635_1[chunk_range])
  my_data4$summidloading[i] <- sum(df2$GCF_002187405_1[chunk_range])
}
#####
ggplot(my_data4,aes(x=SNP_window,y=summaxloading))+geom_line(color="blue1")+
  geom_line(aes(x=SNP_window,y=summinloading),color="firebrick1")+
  geom_line(aes(x=SNP_window,y=summidloading),color="yellow1")+
  ggtitle(paste0("window size ",window_size/1000," kb"))+xlab("SNP window position(Mb)")+ylab("sum(loading)")+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 20))+
  scale_x_continuous(labels = c("0","1","2","3","4","5"))+
  ylim(c(-100,100))

###
ggplot(my_data2,aes(x=log(SNP_distance),y=log(abs(maxPC1strain))))+
  geom_line(color="blue1",linewidth=1)+
  geom_line(aes(x=log(SNP_distance),y=log(abs(midPC1strain))),color="yellow1",linewidth=1)+
  geom_line(aes(x=log(SNP_distance),y=log(abs(minPC1strain))),color="firebrick1",linewidth=1)+
  theme_light()+
  geom_line(data = my_data,aes(x=log(SNP_distance),y=log(abs(Pearson_Correlation_Coefficient))),color="black",linewidth=1)+
  geom_line(data = ldws500,aes(x=log(V2*500),y=log(abs(V3))),color="#C850C0",linewidth=1)+
  ylab("log(Y axis)")+
  xlab("log(SNP distance)")+
  ggtitle("Window Size = 500 bp")+
  annotate("text", x = 6, y = -4, label = "Maximum PC1 strain", color = "blue1", size = 5,hjust = 0, fontface = "bold")+
  annotate("text", x = 6, y = -4.5, label = "Middle PC1 strain", color = "yellow1", size = 5,hjust = 0, fontface = "bold")+
  annotate("text", x = 6, y = -5, label = "Minimum PC1 strain", color = "firebrick1", size = 5,hjust = 0, fontface = "bold")+
  annotate("text", x = 6, y = -5.5, label = "PC1 SNP loading", color = "black", size = 5,hjust = 0, fontface = "bold")+
  annotate("text", x = 6, y = -6, label = "Linkage Disequilibrium", color = "#C850C0", size = 5,hjust = 0, fontface = "bold")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 20))+
  scale_x_continuous(labels = function(x) round(exp(x), 0)) + 
  scale_y_continuous(labels = function(y) round(exp(y), 4)) 

###halfchr
chr1<-read.table("half_matching/halfchrom/chr1pca.eigenvec")
colnames(chr1)[3:22]<-paste0("DNA Fragment 1 PC",c(1:20))
rownames(chr1)<-chr1$V1
chr1$V1<-NULL
chr1$V2<-NULL
chr1val<-read.table("half_matching/halfchrom/chr1pca.eigenval")
colnames(chr1val)<-"val"
rownames(chr1val)<-paste0("DNA Fragment 1 PC", 1:nrow(val))
chr1<-as.data.frame(t(t(chr1)*sqrt(chr1val$val)))
chr1$ID<-rownames(chr1)

chr2<-read.table("half_matching/halfchrom/chr2pca.eigenvec")
colnames(chr2)[3:22]<-paste0("DNA Fragment 2 PC",c(1:20))
rownames(chr2)<-chr2$V1
chr2$V1<-NULL
chr2$V2<-NULL
chr2val<-read.table("half_matching/halfchrom/chr2pca.eigenval")
colnames(chr2val)<-"val"
rownames(chr2val)<-paste0("DNA Fragment 2 PC", 1:nrow(val))
chr2<-as.data.frame(t(t(chr2)*sqrt(chr2val$val)))
chr2$ID<-rownames(chr2)

halfchr<-merge(chr1,chr2)
library(ggpubr)
ggplot(halfchr,aes(x=`DNA Fragment 1 PC1`,y=`DNA Fragment 2 PC1`))+
  geom_point(size=3,alpha=0.5)+
  theme_bw()+
  geom_smooth(method = 'lm', color = 'red',formula = 'y ~ x')+
  stat_cor(data = halfchr,
           size = 6)+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
        axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
        axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）
        axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
        panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1)    # 增粗外框线条
        )

###halfloci
ll<-read.table("half_matching/halfloci/lowlocipca.eigenvec")
colnames(ll)[3:22]<-paste0("Low loading SNPs PC",c(1:20))
rownames(ll)<-ll$V1
ll$V1<-NULL
ll$V2<-NULL
llval<-read.table("half_matching/halfloci/lowlocipca.eigenval")
colnames(llval)<-"val"
rownames(llval)<-paste0("Low loading SNPs PC", 1:nrow(val))
ll<-as.data.frame(t(t(ll)*sqrt(llval$val)))
ll$ID<-rownames(ll)

hl<-read.table("half_matching/halfloci/highlocipca.eigenvec")
colnames(hl)[3:22]<-paste0("High loading SNPs PC",c(1:20))
rownames(hl)<-hl$V1
hl$V1<-NULL
hl$V2<-NULL
hlval<-read.table("half_matching/halfloci/highlocipca.eigenval")
colnames(hlval)<-"val"
rownames(hlval)<-paste0("High loading SNPs PC", 1:nrow(val))
hl<-as.data.frame(t(t(hl)*sqrt(hlval$val)))
hl$ID<-rownames(hl)

halfloci<-merge(hl,ll)
ggplot(halfloci,aes(x=`Low loading SNPs PC1`,y=`High loading SNPs PC1`))+
  geom_point(size=3,alpha=0.5)+
  theme_bw()+
  geom_smooth(method = 'lm', color = 'red',formula = 'y ~ x')+
  stat_cor(data = halfloci,
           size=6)+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）
    axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
    panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1)    # 增粗外框线条
  )

lltree<-read.newick("half_matching/halfloci/lowloci.fasttree")
ggtree(lltree, layout = "unrooted")+
  geom_treescale(offset = 0.01,fontsize=5, linesize=1, x=-0.2, y=-0.2,width=0.01)

hltree<-read.newick("half_matching/halfloci/highloci.fasttree")
ggtree(hltree, layout = "unrooted")+
  geom_treescale(offset = 0.01,fontsize=5, linesize=1, x=0, y=0,width=0.01)

###half samples
geno<-read.table("NR369maf002.matrix_01",check.names = FALSE)
samp1<-read.table("half_matching/halfsample/Half184pca.eigenvec")
colnames(samp1)[3:22]<-paste0("Right-hand PC",c(1:20))
rownames(samp1)<-samp1$V1
samp1geno<-geno[,samp1$V1]
samp1$V1<-NULL
samp1$V2<-NULL
samp1val<-read.table("half_matching/halfsample/Half184pca.eigenval")
colnames(samp1val)<-"val"
rownames(samp1val)<-paste0("Right-hand PC", 1:nrow(samp1val))
samp1<-as.data.frame(t(t(samp1)*sqrt(samp1val$val)))

samp2<-read.table("half_matching/halfsample/Half185pca.eigenvec")
colnames(samp2)[3:22]<-paste0("Left-hand PC",c(1:20))
rownames(samp2)<-samp2$V1
samp2geno<-geno[,samp2$V1]
samp2$V1<-NULL
samp2$V2<-NULL
samp2val<-read.table("half_matching/halfsample/Half185pca.eigenval")
colnames(samp2val)<-"val"
rownames(samp2val)<-paste0("Left-hand PC", 1:nrow(samp2val))
samp2<-as.data.frame(t(t(samp2)*sqrt(samp2val$val)))

samp1_snp <- t(samp1) %*% t(samp1geno)

samp1_samp2<- samp1_snp %*% as.matrix(samp2geno)
samp1_samp2<-as.data.frame(t(samp1_samp2))

halfsample<-cbind(samp1_samp2,samp2)
ggplot(halfsample,aes(x=`Right-hand PC1`,y=`Left-hand PC1`))+
  geom_point(size=3,alpha=0.5)+
  theme_bw()+
  geom_smooth(method = 'lm', color = 'red',formula = 'y ~ x')+
  stat_cor(data = halfsample, label.x = -5000,
           size=6)+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小
    axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
    panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1)    # 增粗外框线条
  )

samp2_snp <- t(samp2) %*% t(samp2geno)

samp2_samp1<- samp2_snp %*% as.matrix(samp1geno)
samp2_samp1<-as.data.frame(t(samp2_samp1))

halfsample2<-cbind(samp2_samp1,samp1)
ggplot(halfsample2,aes(x=`Right-hand PC1`,y=`Left-hand PC1`))+
  geom_point(size=3,alpha=0.5)+
  theme_bw()+
  geom_smooth(method = 'lm', color = 'red',formula = 'y ~ x')+
  stat_cor(data = halfsample2, label.x = -0.2,
           size = 6)+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小
    axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
    panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1)    # 增粗外框线条
  )

###vis half sample
samp1_labeled <- data.frame(ID = row.names(samp1), Source = "Right-hand") 
samp2_labeled <- data.frame(ID = row.names(samp2), Source = "Left-hand")
result_df <- rbind(samp1_labeled, samp2_labeled) 
pca <- read.table("NR369pca.eigenvec")
names(pca)[3:ncol(pca)] <- paste0("PC", c(1:20))
rownames(pca)<-pca$V1
pca$V1<-NULL
pca$V2<-NULL
val<-read.table("NR369pca.eigenval")
colnames(val)<-"val"
rownames(val)<-paste0("PC", 1:nrow(val))
pca<-as.data.frame(t(t(pca)*sqrt(val$val)))
pca$ID<-rownames(pca)

pca_half <-left_join(pca, result_df)

ggplot(pca_half,aes(x=PC1,y=PC2,color=Source,shape=Source))+
  geom_point(size=3)+
  theme_bw()+
  scale_color_aaas()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要 >
    axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
    panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1),    # 增粗外框线条
    legend.title = element_blank(),
    legend.position = c(0.7, 0.85),
    legend.text = element_text(size = 16),
    legend.key.size = unit(2, "lines")
  )

tree.a<-read.newick("NR369_accessions.clean.snp.Fast.tree_nt")
pca_half$label<-pca_half$ID
tree.b<-full_join(tree.a,pca_half,by="label")

ggtree(tree.b,layout = "ape")+
  geom_tippoint(aes(shape = Source,color=Source),size=3)+
  scale_color_aaas()+
  theme(legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.key.size = unit(1, "lines"),
        legend.position = "bottom")+
  geom_treescale(offset = 0.001,fontsize=5, linesize=1, x=0, y=0,width=0.01)
###
###half samples snp loadings
geno<-read.table("NR369maf002.matrix_01",check.names = FALSE)
samp1<-read.table("half_matching/halfsample/Half184pca.eigenvec")
colnames(samp1)[3:22]<-paste0("Right-hand PC",c(1:20))
rownames(samp1)<-samp1$V1
samp1geno<-geno[,samp1$V1]
samp1$V1<-NULL
samp1$V2<-NULL
samp1val<-read.table("half_matching/halfsample/Half184pca.eigenval")
colnames(samp1val)<-"val"
rownames(samp1val)<-paste0("Right-hand PC", 1:nrow(samp1val))
samp1<-as.data.frame(t(t(samp1)*sqrt(samp1val$val)))

samp2<-read.table("half_matching/halfsample/Half185pca.eigenvec")
colnames(samp2)[3:22]<-paste0("Left-hand PC",c(1:20))
rownames(samp2)<-samp2$V1
samp2geno<-geno[,samp2$V1]
samp2$V1<-NULL
samp2$V2<-NULL
samp2val<-read.table("half_matching/halfsample/Half185pca.eigenval")
colnames(samp2val)<-"val"
rownames(samp2val)<-paste0("Left-hand PC", 1:nrow(samp2val))
samp2<-as.data.frame(t(t(samp2)*sqrt(samp2val$val)))

samp1_snp <- t(samp1) %*% t(samp1geno)
samp1_snp <- as.data.frame(t(samp1_snp))
samp2_snp <- t(samp2) %*% t(samp2geno)
samp2_snp <- as.data.frame(t(samp2_snp))

halfloading<-cbind(samp1_snp,samp2_snp)
ggplot(halfloading,aes(x=`Right-hand PC1`,y=`Left-hand PC1`))+
  geom_point(size=3,alpha=0.1)+
  theme_bw()+
  geom_smooth(method = 'lm', color = 'red',formula = 'y ~ x')+
  stat_cor(data = halfloading, label.x = 0,
           size = 6)+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小
    axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
    panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1)    # 增粗外框线条
  )+labs(
    x = "SNP Loading on Right-hand PC1",
    y = "SNP Loading on Left-hand PC1"
  )+xlim(c(-15,20))+ylim(c(-15,20))
halfloading$POS<-as.numeric(rownames(halfloading))
pos2<-load_snp_df[(load_snp_df$PC1>-5)&(load_snp_df$PC1<5),"POS"]
halfloading2<-halfloading[halfloading$POS %in% pos2,]
plot(load_snp_df$POS,load_snp_df$PC1,ylim = c(-5,5))
ggplot(halfloading2,aes(x=`Right-hand PC1`,y=`Left-hand PC1`))+
  geom_point(size=3,alpha=0.1)+
  theme_bw()+
  geom_smooth(method = 'lm', color = 'red',formula = 'y ~ x')+
  stat_cor(data = halfloading2, label.x = 0,
           size = 6)+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小
    axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
    panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1)    # 增粗外框线条
  )+labs(
    x = "SNP Loading on Right-hand PC1",
    y = "SNP Loading on Left-hand PC1",
    title = "c(-5, 5)"
  )
pos3<-load_snp_df[(load_snp_df$PC1>-1)&(load_snp_df$PC1<1),"POS"]
halfloading3<-halfloading[halfloading$POS %in% pos3,]
plot(load_snp_df$POS,load_snp_df$PC1,ylim = c(-1,1))
ggplot(halfloading3,aes(x=`Right-hand PC1`,y=`Left-hand PC1`))+
  geom_point(size=3,alpha=0.1)+
  theme_bw()+
  geom_smooth(method = 'lm', color = 'red',formula = 'y ~ x')+
  stat_cor(data = halfloading3, label.x = 0,
           size = 6)+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小
    axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
    panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1)    # 增粗外框线条
  )+labs(
    x = "SNP Loading on Right-hand PC1",
    y = "SNP Loading on Left-hand PC1",
    title = "c(-1, 1)"
  )
pos4<-load_snp_df[(load_snp_df$PC1>-2)&(load_snp_df$PC1<2),"POS"]
halfloading4<-halfloading[halfloading$POS %in% pos4,]
plot(load_snp_df$POS,load_snp_df$PC1,ylim = c(-2,2))
ggplot(halfloading4,aes(x=`Right-hand PC1`,y=`Left-hand PC1`))+
  geom_point(size=3,alpha=0.1)+
  theme_bw()+
  geom_smooth(method = 'lm', color = 'red',formula = 'y ~ x')+
  stat_cor(data = halfloading4, label.x = 0,
           size = 6)+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小
    axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
    panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1)    # 增粗外框线条
  )+labs(
    x = "SNP Loading on Right-hand PC1",
    y = "SNP Loading on Left-hand PC1",
    title = "c(-2, 2)"
  )

halfloading5<-halfloading[(halfloading$`Right-hand PC1`>-1)&(halfloading$`Right-hand PC1`<1),]
ggplot(halfloading5,aes(x=`Right-hand PC1`,y=`Left-hand PC1`))+
  geom_point(size=3,alpha=0.1)+
  theme_bw()+
  geom_smooth(method = 'lm', color = 'red',formula = 'y ~ x')+
  stat_cor(data = halfloading5, label.x = 0,
           size = 6)+
  theme(
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小
    axis.title.y = element_text(size = 16), # 设置y轴标题的字体大小
    panel.grid = element_blank(),            # 去除所有网格线
    panel.border = element_rect(linewidth = 1)    # 增粗外框线条
  )+labs(
    x = "SNP Loading on Right-hand PC1",
    y = "SNP Loading on Left-hand PC1"
  )

###SNP dists
# 加载必要的库
library(circlize)
library(ComplexHeatmap)

# 读取PCA数据和特征值
pca <- read.table("NR369pca.eigenvec")
eigenval <- scan("NR369pca.eigenval")

# 处理PCA数据，设置行名并移除不需要的列
rownames(pca) <- pca$V1
pca <- pca[,-c(1,2)]
colnames(pca) <- paste0("PC", 1:ncol(pca))

# 读取特征值，并与PCA数据结合
val <- read.table("NR369pca.eigenval")
colnames(val) <- "val"
rownames(val) <- paste0("PC", 1:nrow(val))

# 根据特征值对PCA数据进行缩放
pca <- as.data.frame(t(t(pca) * sqrt(val$val)))

# 设置ID列，并按照PC1排序
pca$ID <- rownames(pca)
pca <- pca[order(pca$PC1),]

# 读取SNP距离数据，并根据PCA的ID排序
snpdist <- read.table("NR369.snp_dists")
snpdist <- snpdist[pca$ID, pca$ID]

# # 定义自定义颜色
# custom_colors <- c("#FF0000", "#800080", "#8A2BE2", "#7FBF7F", "#00FF00")

# 优化颜色映射
color_mapping <- circlize::colorRamp2(
  seq(min(pca$PC1), max(pca$PC1),  length.out = 369), 
  viridis(369, option = "H") 
)

# 创建行注释
row_anno <- rowAnnotation(
  PC1 = pca$PC1, 
  col = list(PC1 = color_mapping),
  show_annotation_name = FALSE,
  show_legend = FALSE  # 不显示图例
)

# 创建列注释
top_anno <- HeatmapAnnotation(
  PC1 = pca$PC1, 
  col = list(PC1 = color_mapping),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    PC1 = list(title = "Strain",
                  at=seq(-0.3,0.4,by=0.2),
                  title_gp = gpar(fontsize = 13),  
                  labels_gp = gpar(fontsize = 12))
  )
)

# 创建热图
heatmap <- Heatmap(as.matrix(snpdist), 
                   name = "SNP Distance",        # 色标名称
                   cluster_rows = FALSE,   # 禁用行聚类
                   cluster_columns = FALSE, # 禁用列聚类
                   show_row_names = FALSE,  # 不显示行名
                   show_column_names = FALSE, # 不显示列名
                   left_annotation = row_anno,   # 添加行注释
                   top_annotation = top_anno ,  # 添加列注释
                   heatmap_legend_param = list(
                     title_gp = gpar(fontsize = 13),  
                     labels_gp = gpar(fontsize = 12) 
                   )
)

# 绘制热图
draw(heatmap, 
     heatmap_legend_side = "right",       # 主图例放右侧
     annotation_legend_side = "right",   # 注释图例放右侧
     merge_legends = TRUE                # 合并图例为单列显示
)

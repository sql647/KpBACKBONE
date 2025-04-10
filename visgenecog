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
geneloading_anno$PC1<-geneloading_anno$PC1^2
geneloading_anno <- geneloading_anno %>%
  mutate(Function = factor(Function, levels = levels(summary_df$Function)))
summary_df <- geneloading_anno %>%
  filter(!is.na(PC1)) %>%
  group_by(COG_category, Function, Category) %>%  # ✅ 加上 Function 和 Category
  summarise(
    mean_PC1sq = mean(PC1, na.rm = TRUE),
    sd_PC1sq   = sd(PC1, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
summary_df <- summary_df %>%
  arrange(Category, desc(mean_PC1sq)) %>%
  mutate(Function = factor(Function, levels = unique(Function)))

ggplot() +
  # 原始点（每个基因的 PC1²）
  geom_jitter(
    data = geneloading_anno,
    aes(x = Function, y = PC1, color = Category),
    width = 0.25, alpha = 0.5, size = 1.5, show.legend = FALSE
  ) +
  # 柱子（平均值）
  geom_col(
    data = summary_df,
    aes(x = Function, y = mean_PC1sq, fill = Category),
    width = 1
  ) +
  # 误差线（SD）
  geom_errorbar(
    data = summary_df,
    aes(
      x = Function,
      ymin = mean_PC1sq - sd_PC1sq,
      ymax = mean_PC1sq + sd_PC1sq
    ),
    width = 0.3, linewidth = 0.6
  ) +
  coord_flip() +
  scale_fill_npg(name = "COG Category") +
  scale_color_npg() +  
  labs(
    x = NULL,
    y = expression("Accessory Genes " * PC1^2),
    title = NULL
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+
  scale_y_continuous(breaks = seq(0, 1000, by = 300))+
  scale_x_discrete(limits = rev)

library(ggplot2)
library(ggsci)

m1<-read.table("hetero1/vcf_pop/pop2maf002.matrix.MajorAMinorT.statRR.avg.matrix.step10.All.stat")
m2<-read.table("hetero2/vcf_pop/pop1maf002.matrix.MajorAMinorT.statRR.avg.matrix.step10.All.stat")

realld<-read.table("../../paper_phd/KP/NR369maf002.matrix.MajorAMinorT.statRR.avg.matrix.step10.All.stat")

m1$group <- "hetero_mutation_rate"
m2$group <- "hetero_HGT_hotspot"

realld$group<-"Real Data"

combined_data <- rbind(m1,m2,realld)

ggplot(data = combined_data) +
  geom_line(aes(x = log(V2 * 10), y = log(V3), color = group), linewidth = 1) +
  xlab("SNP Distance (log scale)") +
  ylab("LD (log scale)") +
  theme_light() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 14),   # 设置x轴刻度标签的字体大小  
    axis.text.y = element_text(size = 14),   # 设置y轴刻度标签的字体大小  
    axis.title.x = element_text(size = 16),  # 设置x轴标题的字体大小（如果需要的话）  
    axis.title.y = element_text(size = 16),  # 设置y轴标题的字体大小（如果需要的话）
    plot.title = element_text(size = 20),
    legend.position =c(0.3,0.25),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20)
  )+
  scale_color_nejm(name = "models",
                  breaks = c("hetero_mutation_rate","hetero_HGT_hotspot","Real Data"))


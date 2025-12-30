####test final state

for (i in c(1:2)) {
  assign(paste0("pop",i),read.table(paste0("./out_vcf/pop",i,"_10000.vcf"), skip = 14, header = TRUE, 
                                    sep = "\t",comment.char = "",check.names = FALSE))
  
}
for (i in c(1:2)) {
  current_df <- get(paste0("pop", i))  # 获取当前的数据框
  current_df <- current_df[,c(2,10:259)]
  colnames(current_df)[2:251] <- paste0("pop",i,"_",c(1:250))
  current_df<-current_df[!duplicated(current_df$POS,fromLast = TRUE),]
  assign(paste0("pop", i), current_df)  # 将修改后的数据框重新赋值
}
df_list <- list(pop1, pop2)
pops <- Reduce(function(x, y) merge(x, y, by = "POS", all = TRUE), df_list)

pops[is.na(pops)] <- 0
pops$`#CHROM`<- 1
pops$ID<-pops$POS
pops$REF<-"A"
pops$ALT<-"T"
pops$QUAL<- 1000
pops$FILTER<-"PASS"
pops$FORMAT<-"GT"
pops$INFO<-"."
pops<-pops[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
              paste0("pop1_",c(1:50)),paste0("pop2_",c(1:50)))]
#pops<-pops[!duplicated(pops$POS,fromLast = TRUE),]

write.table(pops,"vcf_pop/pop2_generation_10000.vcf",
            sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)


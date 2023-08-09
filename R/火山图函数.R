#火山图函数
draw_volcano_plot <- function(exp_data,DEG_list, adj_P_Val = 0.05, hetamap_number=40,logFC_cutoff = 0,  picDir = "火山图") {
 
  if (!file.exists(picDir)) {dir.create(picDir)} 
  library(ggpubr)
  colnames(DEG_list)[1] <- "SYMBOL"
  # 假设你的数据框名称为df
  colnames(DEG_list)[colnames(DEG_list) == "log2FoldChange"] <- "logFC"
  
  # Transform adj.P.Val to -log10(adj.P.Val)
  DEG_list$LogP <- -log10(DEG_list[, grep("adj", names(DEG_list))])#新增-log10p列
  
  # Add Group column based on thresholds
  DEG_list$Group <- "not-significant"
  DEG_list$Group[which((DEG_list[, grep("adj", names(DEG_list))] < adj_P_Val) & (DEG_list$logFC < -logFC_cutoff))] <- "down-regulated"
  DEG_list$Group[which((DEG_list[, grep("adj", names(DEG_list))] < adj_P_Val) & (DEG_list$logFC > logFC_cutoff))] <- "up-regulated"
  # Add top 20 genes as labels
  DEG_list$label <- ""
  top_up_genes <- head(DEG_list$SYMBOL[which(DEG_list$Group == "up-regulated")], 10)
  top_down_genes <- head(DEG_list$SYMBOL[which(DEG_list$Group == "down-regulated")], 10)
  top_genes <- c(as.character(top_up_genes), as.character(top_down_genes))
  DEG_list$label[match(top_genes, DEG_list$SYMBOL)] <- top_genes
  
  # Draw volcano plot
  p <- ggscatter(DEG_list, x = "logFC", y = "LogP", color = "Group", palette = c("#2f5688", "#BBBBBB", "#cc0000"), size = 1, label = DEG_list$label, font.label = 8, repel = T, xlab = "Log2FoldChange", ylab = "-Log10(Adjust P-value)") +
    theme() +
    geom_hline(yintercept = 1.30, linetype = "dashed") +
    if (logFC_cutoff != 0) {
      geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed")
    }
  
  # Save plot to file
  pdf(paste0(picDir, "/火山图.pdf"))
  print(p)
  dev.off()
  
  DEG_list$Group <- factor(DEG_list$Group, levels = c("up-regulated","down-regulated","not-significant"))
  
  
  #########################################        pheatmap      ##################################
  
  library(pheatmap)
  exp_data = as.data.frame(exp_data)
  # 实验组
  exp_data_T <- exp_data %>% dplyr::select(str_which(colnames(.), "_treat")) 
  treatNum <- ncol(exp_data_T) 
  cat("实验组数量:", treatNum, "\n")
  # 正常组
  exp_data_N <- exp_data %>% dplyr::select(str_which(colnames(.), "_con"))
  conNum <- ncol(exp_data_N) 
  cat("正常组数量:", conNum, "\n")
  # 合并数据，正常放前面，肿瘤放后面
  expr <- cbind(exp_data_N, exp_data_T)
  
  top_up_genes <- head(DEG_list$SYMBOL[which(DEG_list$Group == "up-regulated")], hetamap_number/2)
  top_down_genes <- head(DEG_list$SYMBOL[which(DEG_list$Group == "down-regulated")], hetamap_number/2)
  top_genes <- c(as.character(top_up_genes), as.character(top_down_genes))    
  
  annotation_col = data.frame(
  Type = factor(c(rep("con",conNum),rep("treat",treatNum))))
  rownames(annotation_col)=colnames(expr)
  pdf(paste0(picDir,"/前",hetamap_number,"差异基因绘制的热图.pdf"), width = 8, height = 6)
  
  pheatmap(expr[top_genes,], fontsize = 8,
           method="spearman", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
           scale="row", #为基因做scale
           cluster_rows=T,#为基因做聚类
           cluster_cols=F,#为sample做聚类,
           gaps_row=hetamap_number/2,
           gaps_col= conNum,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(20),
           show_colnames=F,show_rownames =T,
           annotation_col = annotation_col,
           # annotation_colors = ann_colors,
           #treeheight_row = "0",treeheight_col = "0",#不画树
           border_color = "NA")
  dev.off()

  
  return(DEG_list)
  
}

ydx_ConsensusCluster <-function(exp_data,genes=t1dm_mci_wgcna_mit,number_of_cluster = 2,pheatmap_width = 10,pheatmap_height =6){
  library(pheatmap)
  library(data.table)
  library(ConsensusClusterPlus)
  library(ComplexHeatmap)
  library(lolR)
  library(impute)
  library(ggplot2)
  library(ggpubr)

  cat("分类算法只针对于疾病组\n")
  exp_data = as.data.frame(exp_data)
  # 实验组
  exp_data_T <- exp_data %>% dplyr::select(grep("_treat", colnames(.)))
  treatNum <- ncol(exp_data_T) 
  cat("实验组数量:", treatNum, "\n")
  
  standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {
    outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
    if (!is.null(halfwidth)) {
      outdata[outdata>halfwidth]=halfwidth
      outdata[outdata<(-halfwidth)]= -halfwidth
    }
    return(outdata)
  }
  # 设置颜色
  red <- "#E53837"
  blue <- "#6AADF0"
  yellow <- "#F7A400"
  heatmap.BLYW <- c("#3F047D","white","#B05900")
  
  
  indata <- exp_data_T[genes,]
  indata <- as.matrix(indata)

  

  cc <- ConsensusClusterPlus(d = indata, 
                             maxK = 10,
                             reps = 1000,
                             pItem = 0.8, 
                             pFeature = 0.8, 
                             clusterAlg = "hc",
                             distance = "pearson", 
                             title = "Consensus Cluster",
                             writeTable = TRUE,
                             plot = "pdf")
  
  
  group <- paste0("Cluster",cc[[number_of_cluster]]$consensusClass); names(group) <- colnames(indata) # 取出三个亚型的结果
  
  
  group <- sort(group) # 按照类排序
  
  plotdata <- standarize.fun(indata[,names(group)],halfwidth = 2) # 绘图数据标准化（z-score并将绝对值超过2的数值截断）
  
  annCol <- data.frame(Cluster = as.character(group), # 构建样本注释
                            row.names = names(group),
                            stringsAsFactors = F)
  annColors <- list("Cluster" = setNames(RColorBrewer::brewer.pal(8, "Set2"), paste0("Cluster", 1:8))) # 构建颜色注释
  
  # 绘制热图
  pdf(paste0("Consensus Cluster", "/聚类热图.pdf"),width = pheatmap_width,height = pheatmap_height)
  plot2= pheatmap(plotdata,
           color = NMF:::ccRamp(x = heatmap.BLYW,n=64), # 原文颜色模版
           annotation_col = annCol[colnames(plotdata),,drop = F],
           annotation_colors = annColors,
           cluster_cols = FALSE, # 样本不聚类，按照亚型顺序排列
           cluster_rows = TRUE, # 行聚类
           border_color = NA,
           treeheight_row = 20, # 修改行树高
           cellheight = 10, # 修改每个单元的高度
           cellwidth = 0.6, # 修改每个单元的宽度
           show_rownames = TRUE, # 显示行名
           show_colnames = FALSE)
  print(plot2)
  dev.off()
  
  
  exp_data_T_rename =t(exp_data_T) 
  merged_df <- cbind(annCol, exp_data_T_rename)
  # 将cluster列的内容添加到现有的行名里面
  rownames(merged_df) <- paste0(rownames(merged_df), "_",merged_df$Cluster)
  # 删除cluster列
  merged_df <- merged_df[, -which(names(merged_df) == "Cluster")]
  
  cluster_name_data = as.data.frame(t(merged_df))
 # colnames(cluster_name_data) <- gsub("_treat", "", colnames(cluster_name_data))
  
  cat("分类算法只针对于疾病组，故里面只有疾病组的信息\n")
  return (list(cluster_table=annCol,cluster_name_data=cluster_name_data ))
  
}
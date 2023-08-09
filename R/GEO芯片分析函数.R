
#基因提取及命名symbol
getGEO_data <- function(GSE_id, GPL_file,symbol_name= "", log_transform = FALSE) {
  options("repos"="https://mirrors.ustc.edu.cn/CRAN/")
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
  current_datetime <- format(Sys.time(), "%Y-%m-%d")
  picDir = paste0(current_datetime,"初步分析" )
  library(readxl)
  library(tidyverse)
  library(GEOquery)
  library(limma) 
  library(affy)
  library(stringr)
  library(FactoMineR)
  library(factoextra)
  library(sva)
  # download and extract GEO data
  gset <- getGEO(GSE_id, destdir = ".", AnnotGPL = TRUE, getGPL = TRUE)
  # process gene expression data
  gse_exp <- exprs(gset[[1]])
  gse_probe_name <- rownames(gse_exp)
  gse_exp_averaged <- avereps(gse_exp)
  # process GPL data
  GPL_data <- Table(getGEO(filename = GPL_file, AnnotGPL = TRUE))
  
  if(symbol_name == ""){
    symbol_name <- grep("gene.*symbol|symbol.*gene", names(GPL_data), ignore.case = TRUE, value = TRUE)
    
  }else{
    symbol_name= symbol_name
  }
  
  print(symbol_name)
  
  
  loc <- match(GPL_data[, 1], gse_probe_name)
  probe_exp <- gse_exp_averaged[loc, ]
  raw_geneid <- as.matrix(GPL_data[, symbol_name])
  index <- which(!is.na(raw_geneid))
  geneid <- raw_geneid[index]
  exp_matrix <- probe_exp[index, ]
  geneid_factor <- factor(geneid)
  gene_exp_matrix <- apply(exp_matrix, 2, function(x) tapply(x, geneid_factor, mean))
  rownames(gene_exp_matrix) <- levels(geneid_factor)
  gse_exp_df <- as.data.frame(gene_exp_matrix)
  gse_exp_df <- na.omit(gse_exp_df)#去除NA列
  # process pData
  gse_pData <- pData(gset[[1]])
  print("数据已转换为symbol名称，并去除NA")

  ex <- gse_exp_df
  qx <- quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  LogC <- (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ||
    any(qx[5:6] > 100)
  if (LogC) {
    ex[ex <= 0] <- NaN
    exprSet <- log2(ex)
    print("需要取log2")
  } else {
    print("无需取log2")
  }
  
  # 取log
  if (log_transform) {
    expr <- log2(expr + 1)
  }
  if (!file.exists(picDir)) {
    dir.create(picDir)
  } 
  p <- boxplot(gse_exp_df, las=1, outline=FALSE,main = "Boxplots")
  p

  pdf(paste0(picDir, "/原始box图看分布.PDF"), height = 9, width = 15)
  boxplot(gse_exp_df, las=1, outline=FALSE,main = "Boxplots")
  dev.off()
  print("已在初步分析文件夹下创建最初的box图，可以查看数据大致分布情况")
  
  gse_exp_df = normalizeBetweenArrays(as.matrix(gse_exp_df))
    pdf(paste0(picDir, "/校正后box图.PDF"), height = 9, width = 15)
    library(RColorBrewer)
    colors =brewer.pal(10,"Paired")
  boxplot(gse_exp_df,col=colors,  las=1, outline=FALSE,main = "Boxplots",notch =TRUE)
  dev.off()
  print("已在初步分析文件夹下创建校正后的box图，可以查看数据大致分布情况")
  
   print("使用了normalizeBetweenArrays函数矫正批次，适用于同一个数据集内矫正")
  
  
  
  print(range(gse_exp_df))
  return(list(pData = gse_pData, exp = gse_exp_df))
}
#删除离群数据
del_data <- function(targets="", expr, outliers) {
  # 删除离群数据
  if (length(targets) == 0){
    # 根据outliers生成正则表达式模式
    pattern <- paste(outliers, collapse = "|")
    # 使用grep函数找到匹配的列名的索引
    cols_to_remove <- grep(pattern, colnames(data), value = TRUE)
    # 选出不包含cols_to_remove变量中列名的所有列
    expr <- data[, !(colnames(data) %in% cols_to_remove)]
  }else{
    # 去除targets中geo_accession列含有离群数据的行
    targets <- targets[!targets$geo_accession %in% outliers, ]
    # 找出targets和expr共同包含的列名
    selected_GSM <- intersect(rownames(targets), colnames(expr))
    # 选出包含selected_GSM变量中列名的所有列
    expr <- dplyr::select(as.data.frame(expr), c(selected_GSM))
  }

  # normalization
  write.table(expr, file = paste0("symbol基因表达矩阵.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
  return(list(targets = targets, expr = expr))
}

#绘制pca
ydx_draw_plot_pca <- function(data, group_var= "") {
  library(ggplot2)
  library(readxl)
  library(tidyverse)
  library(GEOquery)
  library(limma) 
  library(affy)
  library(stringr)
  library(FactoMineR)
  library(factoextra)
  library(sva)
  exp_data=  data
  exp_data = as.data.frame(exp_data)
  # 正常组
  exp_data_N <- exp_data %>% dplyr::select(str_which(colnames(.), "_con"))
  conNum <- ncol(exp_data_N) 
  cat("正常组数量:", conNum, "\n")
  # 实验组
  exp_data_T <- exp_data %>% dplyr::select(str_which(colnames(.), "_treat")) 
  treatNum <- ncol(exp_data_T) 
  cat("实验组数量:", treatNum, "\n")
  # 合并数据，正常放前面，肿瘤放后面
  data <- cbind(exp_data_N, exp_data_T)
 
  group_var=c(rep("con",conNum),rep("treat",treatNum))

  #分组信息

  # 获取时间新建文件夹
  current_datetime <- format(Sys.time(), "%Y-%m-%d")
  picDir = paste0(current_datetime,"PCA")
  if (!file.exists(picDir)) {dir.create(picDir) } 
  
  # 检查是否有缺失值
  missing_vals <- sum(is.na(data))
  cat("Number of missing values: ", missing_vals, "\n")
  
  # 使用正则表达式替换列名
  colnames(data) <- gsub("_con|_treat", "", colnames(data))
  

  # 执行PCA并可视化
  data.pca <- prcomp(data, scale. = TRUE)   
  pcaPredict=predict(data.pca)
  
  
  
  pca <- FactoMineR::PCA(t(data), graph = FALSE)
  
  # gse编号图  
  pca_plot <- factoextra::fviz_pca_ind(pca, 
                           col.ind = group_var,
                           palette = c("#00AFBB", "#E7B800"),
                           addEllipses = TRUE, 
                           legend.title = group_var)

ggsave(paste0(picDir, "/pca_plot编号图.pdf"), pca_plot, width = 10, height = 8)
  # 点图




pca_plot1 <- factoextra::fviz_pca_ind(pca, 
                          geom = "point",
                           col.ind = group_var,
                           palette = c("#00AFBB", "#E7B800"),
                           addEllipses = TRUE, 
                           legend.title = group_var)

ggsave(paste0(picDir, "/pca_plot点图.pdf"), pca_plot1, width = 10, height = 8)
return(pca=pca)
}

#差异计算
ydx_deg <- function(exp_data=counts,is_count =TRUE,show_difgene=TRUE,adjpvalue_cut  = 0.05,logFC_cut = 0) {
  # 将 is_count 设置为逻辑值 TRUE
  library(limma)
  library(dplyr)
  library(tibble)
  library(easyTCGA)
  library(stringr)
  # exp_data=  pre_data$exp_data
  exp_data = as.data.frame(exp_data)
  # 正常组
  exp_data_N <- exp_data %>% dplyr::select(grep("_con", colnames(.)))
  conNum <- ncol(exp_data_N) 
  cat("正常组数量:", conNum, "\n")
  # 实验组
  exp_data_T <- exp_data %>% dplyr::select(grep("_treat", colnames(.)))
  treatNum <- ncol(exp_data_T) 
  cat("实验组数量:", treatNum, "\n")
  # 合并数据，正常放前面，肿瘤放后面
  expr <- cbind(exp_data_N, exp_data_T)
  
  
  #差异分析
  Type=c(rep("con",conNum),rep("treat",treatNum))
  Type=as.factor(Type)
  # 提取差异表达基因列表
  if (is_count) {
    
    expr <- expr[rowMeans(expr)>3,]
    
    DEG_list <- diff_analysis(exprset = expr,
                              group = Type,
                              is_count = TRUE,  # 是 count 数据
                              logFC_cut = logFC_cut,  # 可以直接筛选结果
                              adjpvalue_cut = adjpvalue_cut)
    
    DEG_list_org <- diff_analysis(exprset = expr,
                              group = Type,
                              is_count = TRUE)
    
    cat("已过滤在所有重复样本中小于3的基因，表达量太低也没研究意义")
    #####可视化#####
    edgeR = rownames(DEG_list$deg_edger)
    limma = rownames(DEG_list$deg_limma)
    DESeq2 =as.character( DEG_list$deg_deseq2$genesymbol)
    library(VennDiagram)
    venn.diagram(
      x = list(
        'edgeR' = edgeR,
        'limma' = limma,
        'DESeq2' = DESeq2
      ),
      filename = '三包取差异基因.png',
      col = "black",
      fill = c("dodgerblue", "goldenrod1", "darkorange1"),
      alpha = 0.5,
      cex = 0.8,
      cat.col = 'black',
      cat.cex = 0.8,
      cat.fontface = "bold",
      margin = 0.05,
      main = "三种包的差异表达基因比较",
      main.cex = 1.2
    )
    # 获取交集基因
    co_genes <- intersect(intersect(edgeR, limma), DESeq2)
    co_genes_exp <- counts[co_genes,]
    # 使用append()函数将向量作为子项添加到列表中
    cat("\n添加表达趋势\n")
    ###下面的代码，获得差异表达的基因在三种方法中的表达趋势，并添加到表达矩阵中去
    
    # 将tibble对象转换为数据框
    deg_df <- as.data.frame(DEG_list$deg_deseq2)
    # 使用column_to_rownames函数将"genesymbol"列设置为行名
    rownames(deg_df) <- deg_df$genesymbol
    deg_df = deg_df[co_genes,]#
    deg_limma_1 = DEG_list$deg_limma[co_genes,]
    deg_edger_1 = DEG_list$deg_edger[co_genes,]
    
    co_genes_exp$deg_group <- ifelse(deg_df$log2FoldChange>0,"up","down")
    co_genes_exp$limma_group<- ifelse(deg_limma_1$logFC>0,"up","down")
    co_genes_exp$edgeR_group<- ifelse(deg_edger_1$logFC>0,"up","down")

    DEG_list <- append(DEG_list, list(co_genes))
    DEG_list <- append(DEG_list, list(co_genes_exp))   

    cat("写出共同表达差异的基因矩阵\n")
    write.csv(co_genes_exp,"共同表达差异的基因表达矩阵.csv")

    
    
  } else {
 
    DEG_list <- diff_analysis(exprset = expr,
                              group = Type,
                              is_count = FALSE,  # 不是 count 数据
                              logFC_cut = logFC_cut,  # 可以直接筛选结果
                              adjpvalue_cut = adjpvalue_cut)
    DEG_list_org <- diff_analysis(exprset = expr,
                                  group = Type,
                                  is_count = FALSE)
    cat("DEG_list的基因数目（总的数目）：", nrow(DEG_list), "\n")


  }
  if(show_difgene){
    # 从DEG_list中提取logFC向量
    # 提取DEG_list中的logFC列并进行排序
    sorted_index <- order(-DEG_list$deg_limma$logFC)
    
    # 提取排序后第一个元素对应的symbol
    top_symbol <- DEG_list$deg_limma$genesymbol[sorted_index[1]]
    
    # 提取符合条件的基因的表达量，并添加类型信息
    z <- as.data.frame(t(expr[top_symbol,]))
    z$type <- factor(c(rep("con", conNum), rep("treat", treatNum)))
    # 绘制盒图
    library(ggpubr)
    pdf(paste0("取logFC最高的基因.pdf"))
    plot= ggboxplot(z, x = "type", y = top_symbol,
                    width = 0.6, fill = "type",
                    notch = TRUE, palette = c("#00AFBB", "red","#E7B800"),
                    add = "jitter", shape = "type") + 
      stat_compare_means(aes(group =type))
    print(plot)
    dev.off()
    cat("请注意查看最后一个图形的分组信息，该组是实验组中升高最高基因，确定实验组和对照组是否放反")
  }

 
  
  return(list(DEG_list = DEG_list,DEG_list_org=DEG_list_org))

}
#新的差异表达的矩阵
ydx_select_gene_exp <- function(gene, expr) {
  class(gene)
  if (class(gene) == "data.frame") {
    print("目前正在提前第一列为symbol的基因的表达矩阵")
    genename = gene[,1]
    
    gene_expr <- expr[genename, ]

  } else {
    print("现在输入的是基因列表，正在获取该列表的表达矩阵")
    expr<- as.data.frame(expr)
    gene_expr <- expr[gene, , drop = FALSE]
  }

  # 再次转换为data.frame类型
  gene_expr<- as.data.frame(gene_expr)
  
  write.table(gene_expr, file = paste0("----所需要的基因的表达矩阵.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
  
  
  return(gene_expr)
}

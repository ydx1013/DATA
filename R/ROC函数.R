ydx_draw_roc_curve <- function(genesnames, gene_exp, xunhuan=10, picDir = "ROC_",color = "#d62728") {
  
  library(pROC)
  library(dplyr)
  library(stringr)
  library(limma)
  
  
  data_roc <- as.data.frame(gene_exp)
  # 创建一个名为"ROC"的文件夹用于存储图像
  if (!dir.exists(picDir)) {dir.create(picDir)}
  # 实验组
  exp_data_T <- data_roc %>% dplyr::select(str_which(colnames(.), "_treat")) 
  nT <- ncol(exp_data_T) 
  cat("实验组数量:", nT, "\n")
  
  # 正常组
  exp_data_N <- data_roc %>% dplyr::select(str_which(colnames(.), "_con"))
  nN <- ncol(exp_data_N) 
  cat("正常组数量:", nN, "\n")
  
  # 合并数据，正常放前面，肿瘤放后面
  data_roc <- cbind(exp_data_N, exp_data_T)
  data_roc <- avereps(data_roc)
  # 在picDir中创建一个名为new_folder的文件夹
  
  
if (xunhuan == ""){}else{
  dir.create(file.path(picDir,"1-1循环ROC"))
  for(i in 1:xunhuan){
    cat("目前正在循环第",i,"次\n")
    set.seed(i) # 使用循环变量 i 作为随机种子
    # 从列名中含有con和treat的列中筛选出对照组和实验组的列
    con_cols <- grepl(pattern = "con", colnames(data_roc))
    treat_cols <- grepl(pattern = "treat", colnames(data_roc))
    # 随机选取50%的对照组和实验组样本
    con_sample <- sample(colnames(data_roc)[con_cols], size = length(colnames(data_roc)[con_cols])/2, replace = TRUE)
    treat_sample <- sample(colnames(data_roc)[treat_cols], size = length(colnames(data_roc)[treat_cols])/2, replace = TRUE)
    
    new_df_50 <- data_roc[,c( con_sample, treat_sample)]
    new_cols <- colnames(new_df_50)
    new_df_50 <- data_roc[, -which(colnames(data_roc) %in% new_cols)]
    
    
    # 从原始数据框中根据已选取的对照组和实验组样本，创建新的数据框
    new_df <- data_roc[,c( con_sample, treat_sample)]
    picDir
    
    dir.create(file.path(picDir,"1-1循环ROC",i))
    
    # 创建一个名为"ROC"的文件夹用于存储图像
    if (!dir.exists(picDir)) {dir.create(picDir)}
    
    for (gene in genesnames) {
      data_gene <- t(new_df[gene,,drop=F])
      data_gene <- as.data.frame(data_gene)
      data_gene$group <- ifelse(grepl("treat", rownames(data_gene)), 0, 1)
      y <- data_gene$group
      
      # 绘制ROC曲线
      roc1 <- roc(y, as.numeric(data_gene[, gene]))
      auc_val <- round(as.numeric(gsub("[^0-9.]", "", roc1$auc)), 3)
      ci1 <- ci.auc(roc1, method = "bootstrap")
      ciVec <- as.numeric(ci1)
      
      # 存储ROC曲线图像
      # 在文件夹中创建PDF文件，并将文件存储到指定路径
      pdf(file = file.path(picDir, "1-1循环ROC", i, paste0(gene,"_AUC=", auc_val, "---",  "抽样50.pdf")), width = 5, height = 5)
      plot(roc1, print.auc = TRUE, col = color, legacy.axes = T, main = gene)
      text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f",ciVec[1]), "-", sprintf("%.03f",ciVec[3])), col = color)
      dev.off()
    }
    
    for (gene in genesnames) {
      data_gene <- t(new_df_50[gene,,drop=F])
      data_gene <- as.data.frame(data_gene)
      data_gene$group <- ifelse(grepl("treat", rownames(data_gene)), 0, 1)
      y <- data_gene$group
      
      # 绘制ROC曲线
      roc1 <- roc(y, as.numeric(data_gene[, gene]))
      auc_val <- round(as.numeric(gsub("[^0-9.]", "", roc1$auc)), 3)
      ci1 <- ci.auc(roc1, method = "bootstrap")
      ciVec <- as.numeric(ci1)
      
      # 存储ROC曲线图像
      # 在文件夹中创建PDF文件，并将文件存储到指定路径
      pdf(file = file.path(picDir, "1-1循环ROC", i, paste0(gene,"_AUC=", auc_val, "---",  "剩余50.pdf")), width = 5, height = 5)
      plot(roc1, print.auc = TRUE, col = color, legacy.axes = T, main = gene)
      text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f",ciVec[1]), "-", sprintf("%.03f",ciVec[3])), col = color)
      dev.off()
    }
    
  }
  
  
  
#######################3--7 分##########################
  dir.create(file.path(picDir,"3-7循环ROC"))
  for(i in 1:xunhuan){
    
    set.seed(i) # 使用循环变量 i 作为随机种子
    # 从列名中含有con和treat的列中筛选出对照组和实验组的列
    con_cols <- grepl(pattern = "con", colnames(data_roc))
    treat_cols <- grepl(pattern = "treat", colnames(data_roc))
    # 随机选取30%的对照组和实验组样本
    con_sample <- sample(colnames(data_roc)[con_cols], size = length(colnames(data_roc)[con_cols]) * 0.3, replace = TRUE)
    treat_sample <- sample(colnames(data_roc)[treat_cols], size = length(colnames(data_roc)[treat_cols]) * 0.3, replace = TRUE)
    
    
    # 从原始数据框中根据已选取的对照组和实验组样本，创建新的数据框
    new_df_30 <- data_roc[,c( con_sample, treat_sample)]
    new_cols <- colnames(new_df_30)
    new_df_70 <- data_roc[, -which(colnames(data_roc) %in% new_cols)]
    
    
    
    picDir
    
    dir.create(file.path(picDir,"3-7循环ROC",i))
    
    # 创建一个名为"ROC"的文件夹用于存储图像
    if (!dir.exists(picDir)) {dir.create(picDir)}
    
    for (gene in genesnames) {
      data_gene <- t(new_df_30[gene,,drop=F])
      data_gene <- as.data.frame(data_gene)
      data_gene$group <- ifelse(grepl("treat", rownames(data_gene)), 0, 1)
      y <- data_gene$group
      
      # 绘制ROC曲线
      roc1 <- roc(y, as.numeric(data_gene[, gene]))
      auc_val <- round(as.numeric(gsub("[^0-9.]", "", roc1$auc)), 3)
      ci1 <- ci.auc(roc1, method = "bootstrap")
      ciVec <- as.numeric(ci1)
      
      # 存储ROC曲线图像
      # 在文件夹中创建PDF文件，并将文件存储到指定路径
      pdf(file = file.path(picDir, "3-7循环ROC", i, paste0(gene,"_AUC=", auc_val, "---", "随机30.pdf")), width = 5, height = 5)
      plot(roc1, print.auc = TRUE, col = color, legacy.axes = T, main = gene)
      text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f",ciVec[1]), "-", sprintf("%.03f",ciVec[3])), col = color)
      dev.off()
    }
    for (gene in genesnames) {
      data_gene <- t(new_df_70[gene,,drop=F])
      data_gene <- as.data.frame(data_gene)
      data_gene$group <- ifelse(grepl("treat", rownames(data_gene)), 0, 1)
      y <- data_gene$group
      
      # 绘制ROC曲线
      roc1 <- roc(y, as.numeric(data_gene[, gene]))
      auc_val <- round(as.numeric(gsub("[^0-9.]", "", roc1$auc)), 3)
      ci1 <- ci.auc(roc1, method = "bootstrap")
      ciVec <- as.numeric(ci1)
      
      # 存储ROC曲线图像
      # 在文件夹中创建PDF文件，并将文件存储到指定路径
      pdf(file = file.path(picDir, "3-7循环ROC", i, paste0(gene,"_AUC=", auc_val, "---",  "随机70.pdf")), width = 5, height = 5)
      plot(roc1, print.auc = TRUE, col = color, legacy.axes = T, main = gene)
      text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f",ciVec[1]), "-", sprintf("%.03f",ciVec[3])), col = color)
      dev.off()
    }
  }

  
}


  
  for (gene in genesnames) {
    data_gene <- t(data_roc[gene,,drop=F])
    data_gene <- as.data.frame(data_gene)
    data_gene$group <- ifelse(grepl("treat", rownames(data_gene)), 0, 1)
    y <- data_gene$group
    
    # 绘制ROC曲线
    roc1 <- roc(y, as.numeric(data_gene[, gene]))
    auc_val <- round(as.numeric(gsub("[^0-9.]", "", roc1$auc)), 3)
    ci1 <- ci.auc(roc1, method = "bootstrap")
    ciVec <- as.numeric(ci1)
    
    # 存储ROC曲线图像
    pdf(file = paste0(picDir, "/",  "_AUC=", auc_val,"---",gene, ".pdf"), width = 5, height = 5)
    plot(roc1, print.auc = TRUE, col = color, legacy.axes = T, main = gene)
    text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f",ciVec[1]), "-", sprintf("%.03f",ciVec[3])), col = color)
    dev.off()
  }
}
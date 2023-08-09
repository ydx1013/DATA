ydx_draw_pre_machine_learning <- function(exp_data, targets,group_new="Treat") {
  exp_colnames <- colnames(exp_data)
  matched_rownames <- vector(mode = "character", length = length(exp_colnames))
  targets$group <- tolower(targets$group)
  
  suffix <- ifelse(targets$group %in% c("con", "control", "normal"), "_con", "_treat")
  new_rownames <- paste0(rownames(targets), suffix)
  rownames(targets) <- new_rownames
  
  # 遍历exp_data的列名
  for (i in seq_along(exp_colnames)) {
    # 为当前列名生成可能的行名
    
    con_name <- paste0(exp_colnames[i], "_con")
    treat_name <- paste0(exp_colnames[i], "_treat")
    
    # 在targets的行名中查找匹配的行名
    match_idx <- which(rownames(targets) %in% c(con_name, treat_name))
    
    # 检查是否找到匹配的行名
    if (length(match_idx) > 0) {
      print("匹配成功")
      matched_rownames[i] <- rownames(targets)[match_idx[1]]
    } else {
      print("没有找到匹配的行名，保持原始列名")
      matched_rownames[i] <- exp_colnames[i]
    }
  }
  colnames(exp_data) <- matched_rownames
  
  exp_data <- as.data.frame(exp_data)
  # 实验组
  exp_data_T <- exp_data %>% dplyr::select(str_which(colnames(.), "_treat")) 
  nT <- ncol(exp_data_T) 
  cat("实验组数量:", nT, "\n")
  
  # 正常组
  exp_data_N <- exp_data %>% dplyr::select(str_which(colnames(.), "_con"))
  nN <- ncol(exp_data_N) 
  cat("正常组数量:", nN, "\n")
  
  # 合并数据，CON放前面，Test放后面
  exp_data <- cbind(exp_data_N, exp_data_T)
  exp_data <- avereps(exp_data)

  # 读取输入文件，并将行名中的"-"替换为"_"
  rownames(exp_data) <- sub("-", "_", rownames(exp_data))
  rownames(exp_data) <- sub("/", "_", rownames(exp_data))
  
  
  # 保留首行和列名，修改首行数据
  new_data <- data.frame(t(ifelse(grepl("con", colnames(exp_data)), "Control",
                                  ifelse(grepl("treat", colnames(exp_data)), group_new, 
                                         exp_data[1, ]))))
  colnames(new_data) <- colnames(exp_data)
  new_data =t(new_data)
  # 分组信息
  group <- factor(new_data, levels = c("Control", group_new))
  
  cat("group分组成功，结果如下:", group, "\n")
  cat("数据条目数:", nT+nN, "\n")
  return(list(exp_data = exp_data, group = group))

  
}
# 在函数调用后，将会得到一个新的表格matched_exp_data，它的列名被替换为匹配的行名
ydx_draw_lasso <- function(exp_data, picDir = paste0("机器学习-lasso_", GSE_number)) {
  
  # 设置随机种子以保证结果的可重复性
  set.seed(123)
  
  # 引入所需的包
  library(glmnet)


  
  # 将数据转置，以便样本在列上
  rt = t(exp_data)
  
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
  
  # 构建LASSO模型
  x = as.matrix(rt)
  y = sub(".*?(con|treat).*", "\\1", grep("(con|treat)", row.names(rt), value = TRUE))
  fit = glmnet(x, y, family="binomial", alpha=1)
  
  # 绘制lambda值和交叉验证误差之间的曲线图，并保存为pdf文件
  pdf(paste0(picDir, "/lambda.pdf"))
  plot(fit, xvar = "lambda", label = TRUE)
  dev.off()
  
  # 进行10倍交叉验证，并生成交叉验证误差图
  cvfit = cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
  pdf(paste0(picDir, "/cvfit.pdf"),width=6,height=5.5)
  plot(cvfit)
  dev.off()
  
  # 提取拟合LASSO模型后选定的特征基因，并保存到txt文件中
  coef = coef(fit, s = cvfit$lambda.min)
  index = which(coef != 0)
  lassoGene = row.names(coef)[index]
  lassoGene = lassoGene[-1] # 删除响应变量
  write.table(lassoGene, file=paste0(picDir, "/LASSO.gene.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  
  # 提取筛选的特征基因的表达数据，并将其保存到txt文件中
  #输出重要基因的表达量
  sigExp <- t(rt[,lassoGene])
  lassoexp <- rbind(ID = colnames(sigExp), sigExp)
  
  
  write.table(lassoexp, file = paste0(picDir, "/LASSO.geneExp.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  
  
  
  
  return(list(lassoexp=lassoexp, lassoGene=lassoGene))
}
ydx_draw_Decision_Tree <- function(data,importantvlue=2, picDir = paste0("机器学习-决策树_", GSE_number)) {
  library(readr); library(VIM); library(caret); library(rpart); library(rpart.plot); library(Metrics); library(stringr); library(tibble); library(bitops); library(rattle); library(RColorBrewer); library(tidyverse); library(limma); library(pheatmap); library(visNetwork); library(ggpol); library(ggplot2); library(sparkline)
  # 转置数据
  data = t(data)
  
  # 如果目录不存在，则创建
  if (!file.exists(picDir)) {
    dir.create(picDir)
  } 
  
  # 获取组列
  group <- sub(".*?(con|treat).*", "\\1", grep("(con|treat)", row.names(data), value = TRUE))
  
  # 显示数据维度和列名
  cat("数据维度：", dim(data), "\n")
  cat("列名：", colnames(data), "\n")
  
  # 执行数据聚合
  aggr(data)
  
  # 指定随机数种子以保证可重复性
  set.seed(123)
  
  # 将数据转换为data.frame格式
  data2 <- as.data.frame(data)
  
  # 建立决策树模型
  mod1<-rpart(as.factor(group)~.,data = data2,method = "class")
  # 可视化决策树
  pdf(paste0(picDir, "/decision_tree前.pdf"))
  rpart.plot(mod1)
  dev.off()
  
  # 显示变量重要性
  importances <- varImp(mod1)
  jcsimportances=as.matrix(importances %>%
                             arrange(desc(Overall)))
  
  
  JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
  JCSGenes=names(JCSGenes[JCSGenes>importantvlue])
  # 显示模型的最低CP值
  cat("模型最低CP值：", mod1$cp, "\n")
  
  # 绘制复杂度参数图
  pdf(paste0(picDir, "/plotcp.pdf"))
  plotcp(mod1)
  dev.off()
  
  # 用最低CP值构建优化后的决策树模型
  mod1_optimized <- rpart(as.factor(group) ~ ., data = data2, method = "class", cp = 0.00028)
  
  # 将可视化的决策树保存到文件
  library(visNetwork)
  
  
  visTree(mod1_optimized, main="Decision Tree", height="600px",
          colorY=c("greenYellow","hotPink","yellow"), legendWidth=0.2, legendNcol=2,
          nodesFontSize=16, edgesFontSize=10)
  # 可视化决策树
  pdf(paste0(picDir, "/decision_tree后.pdf"))
  rpart.plot(mod1_optimized)
  dev.off()
  write.table(JCSGenes, file=paste0(picDir, "/Decision_Tree.gene_names.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  

  
  #输出重要基因的表达量
  sigExp <- t(data[,JCSGenes])
  JCSGenesexp <- rbind(ID = colnames(sigExp), sigExp)
  
  
  write.table(JCSGenesexp, file = paste0(picDir, "/Decision_Tree.gene_geneExp.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  # 返回优化后的决策树模型
  
  return(list(JCSGenes = JCSGenes, JCSGenesexp = JCSGenesexp))
}
ydx_draw_randomForest <- function(exp_data, importantvlue=2,picDir = paste0("机器学习-randomForest_", GSE_number)) {
  #引用包
  library(randomForest)
  set.seed(123456)
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
  #读取输入文件
  data <- t(exp_data)
  group <- sub(".*?(con|treat).*", "\\1", grep("(con|treat)", row.names(data), value = TRUE))
  
  #随机森林树
  rf <- randomForest(as.factor(group)~., data=data, ntree=1000)
  pdf(file = paste0(picDir, "/森林.pdf"), width=6, height=6)
  plot(rf, main="Random forest", lwd=2)
  dev.off()
  
  #找出误差最小的点
  optionTrees <- which.min(rf$err.rate[,1])
  optionTrees
  rf2 <- randomForest(as.factor(group)~., data=data, ntree=optionTrees)
  
  #查看基因的重要性
  importance <- importance(x=rf2)
  
  #绘制基因的重要性图
  pdf(file = paste0(picDir, "/GeneIm.pdf"), width=6.2, height=5.8)
  varImpPlot(rf2, main="")
  dev.off()
  
  #挑选疾病特征基因
  rfGenes <- importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
  rfGenes2 <- names(rfGenes[rfGenes>0])      #挑选重要性评分大于0的基因
  rfGenes <- names(rfGenes[rfGenes>importantvlue])     #挑选重要性评分大于importantvlue的基因

  write.table(rfGenes, file = paste0(picDir, "/随机森林Genes.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(rfGenes2, file = paste0(picDir, "/所有评分大于0的—Genes.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  #输出重要基因的表达量
  sigExp <- t(data[,rfGenes])
  sigExpOut <- rbind(ID = colnames(sigExp), sigExp)
  write.table(sigExpOut, file = paste0(picDir, "/imGeneExp.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  return(list(genename=rfGenes, Geneexp=sigExpOut))
  
  }
ydx_draw_GBM <- function(exp_data,importantvlue=2, picDir = paste0("机器学习-GBM_", GSE_number)) {
  packages <- c("readr", "VIM", "caret", "rpart", "rpart.plot", "Metrics", "stringr", "tibble", "bitops", "rattle", "RColorBrewer", "tidyverse", "limma", "pheatmap", "visNetwork", "ggpol", "ggplot2", "sparkline", "randomForest", "venn", "dplyr", "caret", "glmnet", "xgboost", "DALEX", "gbm", "VennDiagram", "neuralnet", "NeuralNetTools")
  lapply(packages, require, character.only = TRUE)
  data <- t(exp_data)
  group <- sub(".*?(con|treat).*", "\\1", grep("(con|treat)", row.names(data), value = TRUE))
  
  set.seed(1234)
  metric <- "RMSE"
  myControl <- trainControl(method = "cv", number = 5)
  
  # Fitting model
  fitControl <- trainControl(method = "repeatedcv", number = 4, repeats = 4)
  fit <- train(x = data, y = as.factor(group), method = "gbm", trControl = fitControl, verbose = FALSE)
  
  # 绘制基因重要性梯度图
  importances <- varImp(fit)
  importance <- as.data.frame(importances$importance)
  #绘制基因重要性梯度图



  GBMimportance <- as.matrix(importances$importance)
  GBMGenes=GBMimportance[order(GBMimportance[,"Overall"], decreasing = TRUE),]
  GBMGenes=names(GBMGenes[GBMGenes>importantvlue])
  
  
  
  
  
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
  # 输入完上面那串代码后，显示的结果就是GBM结果，将他们复制出来。
  # 重要性筛选区域设置多少可以自己定，我这里是只要重要性不为0都可以
  # 删除为重要性为0的gene后重新导入
  a <- importance
  varimpdf <- data.frame(var = row.names(a), impor = a[,1])
  
  # 绘制图形并保存为 PDF
  ggplot(varimpdf,aes(x = reorder(var,-impor), y = impor))+
    geom_col(colour = "lightblue",fill = "lightblue")+
    labs(title="Feature gene importance (Gradient Boosting Machine)", x="",y = "importance")+
    theme(plot.title = element_text(size=12,hjust=0.5))+
    theme(axis.text.x = element_text(size = 5))+
    theme(axis.text.y = element_text(size = 12))+
    theme(axis.text.x = element_text(angle = 50,vjust = 0.85,hjust = 0.75))
  ggsave(paste0(picDir, "/GBM.pdf"), width = 15, height = 8) #保存宽度15英寸和高度8英寸,

  #输出重要基因的表达量
  sigExp <- t(data[,GBMGenes])
  GBMGenesexp <- rbind(ID = colnames(sigExp), sigExp)
  
  write.table(GBMGenes, file=paste0(picDir, "/GBMGenes_names.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  
  write.table(GBMGenesexp, file = paste0(picDir, "/GBMGenes_geneExp.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  # 返回优化后的决策树模型
  
  return(list(GBMGenes = GBMGenes, GBMGenesexp = GBMGenesexp))
}
ydx_draw_XGboost <- function(exp_data, importantvlue=0,method = "xgbTree",picDir = paste0("机器学习-XGboost_", GSE_number)) {
  packages <- c("readr", "VIM", "caret", "rpart", "rpart.plot", "Metrics", "stringr", "tibble", "bitops", "rattle", "RColorBrewer", "tidyverse", "limma", "pheatmap", "visNetwork", "ggpol", "ggplot2", "sparkline", "randomForest", "venn", "dplyr", "caret", "glmnet", "xgboost", "DALEX", "gbm", "VennDiagram", "neuralnet", "NeuralNetTools")
  lapply(packages, require, character.only = TRUE)
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
  data=t(exp_data)
  group <- sub(".*?(con|treat).*", "\\1", grep("(con|treat)", row.names(data), value = TRUE))
  # Fitting model(用caret实现)
  TrainControl <- trainControl( method = "repeatedcv", number = 10, repeats = 4)
  model<- train(x=data,y=as.factor(group),  method = method, trControl = TrainControl,verbose = FALSE)
  plot(varImp(model))
  # 把基因重要性绘制成图片
  
  pdf(paste0(picDir, "/重要性基因排序.pdf"))
  plot(varImp(model))
  dev.off()
  

  importance <- varImp(model)
  XGBimportant <- as.matrix(importance$importance)
  XGBGenes=XGBimportant[order(XGBimportant[,"Overall"], decreasing = TRUE),]
  XGBGenes=names(XGBGenes[XGBGenes>importantvlue])
  
  head(importance)
  
  important <- as.data.frame(importance$importance) 
  a<-important
  varimpdf <- data.frame(var = row.names(a),
                         impor = a[,1])
  
  ggplot(varimpdf,aes(x = reorder(var,-impor), y = impor))+
    geom_col(colour = "lightblue",fill = "lightblue")+
    labs(title="Feature gene importance (XGBoost)", x="",y = "importance")+
    theme(plot.title = element_text(size=12,hjust=0.5))+
    theme(axis.text.x = element_text(size = 3))+
    theme(axis.text.y = element_text(size = 12))+
    theme(axis.text.x = element_text(angle = 50,vjust = 0.85,hjust = 0.75))
  pdf <- "/XGboost.pdf" #PDF文件名 可以包含任何后缀
  ggsave(paste0(picDir, pdf), width = 15, height = 8) #注意设置宽度和高度,以PDF推荐的规格保存
  
  
  
  write.table(varimpdf, file=paste0(picDir, "/XGboost.varimpdf.gene_names.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  
  write.table(XGBGenes, file=paste0(picDir, "/XGboost.gene_names.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  
  
  #输出重要基因的表达量
  sigExp <- t(data[,XGBGenes])
  XGBexp <- rbind(ID = colnames(sigExp), sigExp)
  
  
  write.table(XGBexp, file = paste0(picDir, "/XGboost.gene_geneExp.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  # 返回优化后的决策树模型
  
  return(list(XGBGenes = XGBGenes, XGBexp = XGBexp))
  
  
  
  
  
  
  
  
}
ydx_draw_SVM_RFE <- function(exp_data, picDir = paste0("机器学习SVM_", GSE_number)) {
  
  library(e1071)
  library(kernlab)
  library(caret)
  
  set.seed(123)
  
  data <- t(exp_data)
  group <- sub(".*?(con|treat).*", "\\1", grep("(con|treat)", row.names(data), value = TRUE))
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
  #SVM-RFE分析
  Profile <- rfe(x = data,
                 y = as.numeric(as.factor(group)),
                 sizes = c(2, 4, 6, 8, seq(10, 40, by = 3)),
                 rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
                 methods = "svmRadial")
  
  #绘制图形
  pdf(file = file.path(picDir, "SVM-RFE.pdf"), width = 6, height = 5.5)
  par(las = 1)
  x <- Profile$results$Variables
  y <- Profile$results$RMSE
  plot(x, y, xlab = "Variables", ylab = "RMSE (Cross-Validation)", col = "darkgreen")
  lines(x, y, col = "darkgreen")
  #标注交叉验证误差最小的点
  wmin <- which.min(y)
  wmin.x <- x[wmin]
  wmin.y <- y[wmin]
  points(wmin.x, wmin.y, col = "blue", pch = 16)
  text(wmin.x, wmin.y, paste0('N=', wmin.x), pos = 2, col = 2)
  dev.off()
  
  #输出选择的基因
  featureGenes <- Profile$optVariables
  write.table(featureGenes, file = file.path(picDir, "SVM-RFE.gene.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  #输出重要基因的表达量
  sigExp <- t(data[,featureGenes])
  SVM_exp <- rbind(ID = colnames(sigExp), sigExp)
  write.table(SVM_exp, file = paste0(picDir, "/SVM_exp_geneExp.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  # 返回优化后的决策树模型
  
  return(list(SVMGenes = featureGenes, SVM_exp = SVM_exp))
}






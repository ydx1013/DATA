ydx_draw_Nomogram <- function(gene, expr_data,color="blue") {
  # 加载所需库
  library(tidyverse)
  library(rms)
  
  # 检查输入的基因数据类型
  if (class(gene) == "data.frame") {
    # 如果输入为基因表达矩阵，提取第一列作为基因名称
    genename <- gene[, 1]
    gene_expr <- expr_data[genename, ]
  } else {
    # 如果输入为基因列表，获取对应的表达矩阵
    expr_data <- as.data.frame(expr_data)
    gene_expr <- expr_data[gene, , drop = FALSE]
  }
  
  # 转置基因表达矩阵并添加状态列
  gene_expr <- t(gene_expr)
  gene_expr<- as.data.frame(gene_expr)
  gene_expr$status <- ifelse(grepl("_con", rownames(gene_expr)), 0, 1)
  rt <- as.data.frame(gene_expr)
  
  # 打包数据
  ddist <- datadist(rt)
  options(datadist = "ddist")
   
  # 构建逻辑回归模型
  formula <- as.formula(paste("status ~", paste(gene, collapse = "+")))
  fit <- lrm(formula, data = rt, x = TRUE, y = TRUE)
  
  # 构建Nomogram
  nom <- nomogram(fit, fun = plogis, fun.at = c(0.001, 0.01, 0.2, 0.5, 0.8, 0.95, 0.999),
                  lp = FALSE, funlabel = "Risk of IDD")
  
  # 绘制并保存Nomogram为PDF
  pdf("Nomogram.pdf",width = 8,height = 6)
  plot(nom)
  dev.off()
  
  # 计算C-index
  Cindex <- rcorrcens(status ~ predict(fit), data = rt)
  
  # 绘制ROC曲线并保存为PDF
  pdf("ROC_curve.pdf",width = 6,height = 6)
  library(pROC)
  gfit <- roc(status ~ predict(fit), data = rt)
  auc_val <- round(as.numeric(gsub("[^0-9.]", "", gfit$auc)), 3)
  ci1 <- ci.auc(gfit, method = "bootstrap")
  ciVec <- as.numeric(ci1)
  plot(gfit,
       print.auc = TRUE,
       main = "ROC CURVE",
       col = color,
       print.thres.col = "black",
       identity.col = color,
       identity.lty = 1,
       identity.lwd = 1,
       legacy.axes = TRUE)
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f",ciVec[1]), "-", sprintf("%.03f",ciVec[3])), col = color)
  dev.off()
  
  # 绘制Calibration plot并保存为PDF
  pdf("Calibration_plot.pdf",width = 6,height = 6)
  cal <- calibrate(fit, method = "boot", B = 1000)
  plot(cal,
       xlab = "Nomogram-predicted probability of nonadherence",
       ylab = "Actual diagnosed nonadherence (proportion)",
       sub = FALSE)
  dev.off()
  
  # 输出C-index
  Cindex
  
  # 返回结果
  return(list(nomogram = nom, cindex = Cindex))
}


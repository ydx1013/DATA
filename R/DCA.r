dca <- function(expr = pre_data$exp_data, gene_names=c("UQCRH","LYPLAL1","ACADVL","HSPE1")){
  library(rmda)
  
  
  surv <- function(datExpr){
    traitData <- matrix(0, nrow = ncol(datExpr), ncol = 1,
                        dimnames = list(colnames(datExpr), c("outcome")))
    # # 为新矩阵填写值
    traitData[grepl("treat", rownames(traitData)), "outcome"] <- 1
    traitData[grepl("con", rownames(traitData)), "outcome"] <- 0
    traitData <- as.data.frame(traitData)
  }
  
  

  print("现在输入的是基因列表，正在获取该列表的表达矩阵")
  expr<- as.data.frame(expr)
  gene_expr <- expr[gene_names, , drop = FALSE]
  gene_expr= t(gene_expr)
  dcaData1 = surv(expr)
  
  data = cbind(dcaData1,gene_expr)
  
 
  Full_model <- decision_curve(outcome ~ ACADVL+LYPLAL1+HSPE1+UQCRH, # R语言里常见的公式类型
                           data = data, 
                           study.design = "cohort", # 选择研究类型
                           bootstraps = 50 # 重抽样次数
    )
  ACADVL <- decision_curve(outcome ~ ACADVL, # R语言里常见的公式类型
                          data = data, 
                          study.design = "cohort", # 选择研究类型
                          bootstraps = 50 # 重抽样次数
  )
  LYPLAL1 <- decision_curve(outcome ~ LYPLAL1, # R语言里常见的公式类型
                             data = data, 
                             study.design = "cohort", # 选择研究类型
                             bootstraps = 50 # 重抽样次数
  )
  HSPE1 <- decision_curve(outcome ~ HSPE1, # R语言里常见的公式类型
                             data = data, 
                             study.design = "cohort", # 选择研究类型
                             bootstraps = 50 # 重抽样次数
  )
  UQCRH <- decision_curve(outcome ~ UQCRH, # R语言里常见的公式类型
                             data = data, 
                             study.design = "cohort", # 选择研究类型
                             bootstraps = 50 # 重抽样次数
  )
  # 画图
  plot_decision_curve(list(Full_model,ACADVL,LYPLAL1,HSPE1,UQCRH), curve.names = c("Full_model","ACADVL","LYPLAL1","HSPE1","UQCRH"),
                      cost.benefit.axis = F, # 是否需要损失：获益比 轴
                      confidence.intervals = "none" # 不画可信区间
  )
  
  summary(UQCRH)
  
  
}
 

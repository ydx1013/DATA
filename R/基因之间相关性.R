# 绘制单个基因与另一个基因的相关性图
# 参数：
# gene1: 第一个基因的名称
# gene2: 第二个基因的名称
# exp: 基因表达矩阵，每行代表一个基因，每列代表一个样本
ydx_draw_single_gene_relation <- function(gene1="TSPAN6", gene2="TNMD", exp,color = "red") {
  # 提取gene1和gene2在各个样本中的表达值
  x1 <- as.numeric(exp[gene1,])
  y1 <- as.numeric(exp[gene2,])
  
  # 计算两个基因的相关性系数
  corr1 <- cor(x1, y1)
  corr1 <- round(corr1, 5)
  
  # 绘制相关性图，并加入回归线
  pdf(paste0(gene1,"---",gene2,"---","Singlecor.pdf"))
  plot(x1, y1, 
       type="p", # 点的类型,"p" for points
       pch=12, # 点的形状
       cex=.8, 
       col=color, # 点的颜色
       main=paste("Cor r value = ", corr1), # 标题
       xlab=paste(gene1, "expression"), # X轴标签
       ylab=paste(gene2, "expression")) # Y轴标签
  abline(lm(y1 ~ x1, data = as.data.frame(exp))) # 绘制回归线
  dev.off()
}

ydx_draw_gene_to_genes_relation <- function(gene1="ALKBH5", genes, exp, method = "pearson") {
  library(ggstatsplot)
  library(ggside)
  library(ggsci)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tibble)
  
  exp=t(Muldat)
  
  
  exp = exp[genes,]
  # 将基因表达矩阵转换为数据框类型（data.frame）
  Muldat <- as.data.frame(t(exp))

  # 提取gene1在各个样本中的表达值
  y <- as.numeric(Muldat[, gene1])
  
  # 创建一个空的数据框，用于存储相关性系数和P值
  cor_data <- data.frame(genes)
  i = 1
  for (genes_i in genes) { # 对于每个基因genes_i，计算其与gene1之间的相关性系数
    test <- cor.test(as.numeric(Muldat[, genes_i]),# x = Muldat的每一列
                     y, # y已经固定好了，就是gene1那一列
                     method =method, # 统计方法
                     conf.level = 0.95 # 置信区间
    )
    # 将相关性系数和P值存储到cor_data数据框中
    cor_data[i, 2] <- test$estimate # test$estimate 赋值给cor_data的第二列
    cor_data[i, 3] <- test$p.value # test$p.value 赋值给cor_data的第三列
    i = i+1
  }
  # 设置数据框的列名
  names(cor_data) <- c("symbol", "correlation", "pvalue")

  
  pdf("MultiCor1.pdf")
  ggscatterstats(data = Muldat, #要分析的数据
                 y = METTL3, #设置Y轴
                 x = WTAP,#设置X轴
                 margins = "both",#是否显示 边缘，默认为true                                      
                 xfill = "red", #x轴边缘图形的颜色
                 yfill = "blue", #y轴边缘图形的颜色
                 marginal.type = "violin",#在图片坐标轴边缘添加图形类型
                 title = "Relationship between FTO and ZC3H13")#标题
  dev.off()
  
  
  
  
  gene1="ALKBH5"
  genes= c("ALKBH5","METTL3","HNRNPC","YTHDF2","WTAP","KIAA1429")
  # 循环绘制散点图并保存到pdf文件中
  for (genes_2 in genes1) {


    # 提取gene1和gene2在各个样本中的表达值
    pdf(paste0(gene1, "---", genes_i, "---", "Multi Cor1.pdf"))
    ggscatterstats(data = Muldat, # 要分析的数据
                   y = Muldat[[as.name(gene1)]], # 设置Y轴
                   x = Muldat[[as.name(genes_2)]], # 设置X轴
                   margins = "both", # 是否显示边缘，默认为true
                   xfill = "red", # x轴边缘图形的颜色
                   yfill = "blue", # y轴边缘图形的颜色
                   marginal.type = "violin",
                   title = paste0("Relationship between ", gene1, " and ", genes_i)
    )
    dev.off()
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
}

ydx_draw_mult_generelation <- function(genes,exp,main_name= "",picDir = paste0("gene-relation_", GSE_number)){
  exp = exp[genes,]
  exp =t(exp)
  cell=exp
  # 创建一个名为"ROC"的文件夹用于存储图像
  if (!dir.exists(picDir)) {
    dir.create(picDir)
  }
  library(corrplot)
  #使用cor函数和cor.mtest进行相关性分析
  corr <- cor(exp)#计算相关系数
  cor_res <- cor.mtest(exp,#计算p值
                       conf.level = 0.95)#置信区间
  
  
  pdf(file = paste0(picDir, "/","gene-corheatmap.pdf"))
  corrplot(corr,
           method = "pie",#相关性矩阵展示的图形
           order = "AOE",#相关性矩阵图的排序方法
           addCoef.col = "black",#为相关系数添加颜色
           tl.col="black",#设置文本标签的颜色
           p.mat = cor_res$p,#设置P值
           sig.level = 0.05,#设置p值的显著水平，小于0.001为显著p值
           insig = "pch",#对非显著p值打叉，也就是因为p值过大而舍去
           number.cex = 1,#写入相关系数的文本参数
           type = "lower")#右上角显示相关性图
  dev.off()
  
  
  
  library(PerformanceAnalytics)
  pdf(file = paste0(picDir, "/","基因与基因之间相关性图-1.pdf"))
  chart.Correlation(exp,
                    histogram = TRUE,#对角线是否展示histogram
                    pch=19)#散点图打点样式
  dev.off()
  pdf(file = paste0(picDir, "/","基因与基因之间相关性图-2.pdf"))
  chart.Correlation(exp,
                    histogram = F,#对角线是否展示histogram
                    pch=19)#散点图打点样式
  dev.off()


  library(corrgram)
  pdf(file = paste0(picDir, "/","基因与基因之间相关性图-3.pdf"))
  corrgram(cell, 
           order=TRUE, 
           lower.panel = panel.shade, #对角线下方
           upper.panel = panel.pie,#对角线上方，饼图
           text.panel = panel.txt,#设置对角线显示的内容
           cor.method = "pearson",
           main = main_name)
  dev.off()
  
  pdf(file = paste0(picDir, "/","基因与基因之间相关性图-4.pdf"))
  corrgram(cell, 
           order=TRUE, 
           lower.panel = panel.ellipse, #变化了参数，椭圆
           upper.panel = panel.conf,#变化了参数，相关系数
           text.panel = panel.txt,
           cor.method = "pearson",
           main = main_name)
  dev.off()
  pdf(file = paste0(picDir, "/","基因与基因之间相关性图-5.pdf"))
  corrgram(cell, 
           order=TRUE, 
           lower.panel = panel.pie, #饼图
           upper.panel = panel.pts,# 类似于散点图
           text.panel = panel.txt,
           cor.method = "spearman",
           main = main_name)
  dev.off()
  
  #熟悉的corrplot包画相关系数图
  library(corrplot)
  cell = cor(cell)#需要先计算相关系数
  pdf(file = paste0(picDir, "/","基因与基因之间相关性图-6.pdf"))
  corrplot(cell, 
           method = "circle", 
           order = "AOE", 
           tl.pos = "d")
  dev.off()
  #ggcorrplot包画相关系数图
library(ggcorrplot)
cell <- round(cor(cell), 2)
fix(cell)
pdf(file = paste0(picDir, "/","基因与基因之间相关性图-7.pdf"))  
ggcorrplot(cell, 
           hc.order = TRUE, 
           outline.color = "white")
dev.off()

pdf(file = paste0(picDir, "/","基因与基因之间相关性图-8.pdf"))
ggcorrplot(cell, 
           hc.order = TRUE, #使用hc.order进行排序
           type = "lower", #图片位置
           outline.color = "white",#轮廓颜色
           lab = TRUE,#true为在图上添加相关系数
           ggtheme = ggplot2::theme_gray, #指ggplot2函数对象，默认值为thememinimal
           colors = c("#6D9EC1", "white", "#E46726"))#设置颜色
dev.off()

library(PerformanceAnalytics)
pdf(file = paste0(picDir, "/","基因与基因之间相关性图-9.pdf"))
chart.Correlation(cell,
                  histogram = F,
                  pch=1)
dev.off()


#ggcorrplot包画相关系数图
library(ggcorrplot)
cell <- round(cor(cell), 2)

pdf(file = paste0(picDir, "/","基因与基因之间相关性图-10.pdf"))
ggcorrplot(cell, 
           hc.order = TRUE, 
           outline.color = "white")
dev.off()

pdf(file = paste0(picDir, "/","基因与基因之间相关性图-11.pdf"))
ggcorrplot(cell, 
           hc.order = TRUE, #使用hc.order进行排序
           type = "lower", #图片位置
           outline.color = "white",#轮廓颜色
           lab = TRUE,#true为在图上添加相关系数
           ggtheme = ggplot2::theme_gray, #指ggplot2函数对象，默认值为thememinimal
           colors = c("#6D9EC1", "white", "#E46726"))#设置颜色
dev.off()

 
  
  
}
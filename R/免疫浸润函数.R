# Cibersort：(非肿瘤可用)   是一种常用的免疫浸润分析方法，该方法基于已知参考数据集，默认提供22种免疫细胞亚型的基因表达特征集：LM22。其通过不同免疫细胞中标志基因的差异表达分析出样品中各种免疫细胞的种类和分布，可以用于研究不同样品的免疫细胞种类的差异情况。
# 
# xCell：(非肿瘤可用) 是基于 ssGSEA 的方法，可根据 64 种免疫细胞和基质细胞类型的基因表达数据进行细胞类型富集分析。该方法可以估计64中免疫细胞类型的丰度分数，包括适应性和先天免疫细胞、造血祖细胞、上皮细胞、细胞外基质细胞。
# xCell 是基于 ssGSEA 的方法，可根据 64 种免疫细胞和基质细胞类型的基因表达数据进行细胞类型富集分析。xCell 使用表达水平排名而不是实际值，因此归一化没有影响，但是输入数据需要 normalization 格式（RPKM/FPKM/TPM/RSEM）。
#
# IPS：通过免疫表型评分（immunophenoscore, IPS）分别计算四种不同免疫表型的评分（抗原呈递、效应细胞、抑制性细胞、检查点），IPS z-score为四者的整合，且IPS z-score越高，样本免疫原性越强。
#免疫表型评分(IPS)可用于预测对免疫检查点抑制剂的反应

# ESITMATE：(肿瘤可用)研究免疫浸润的一种常用算法，可以通过基因表达量为肿瘤基质细胞和免疫细胞浸润结果评分，可以用于研究不同样品间免疫浸润的情况。利用癌症样本转录谱来推断肿瘤细胞的含量，以及免疫细胞和基质细胞的浸润程度。其他方法不同的是：（1）除了免疫细胞，还能分析肿瘤细胞纯度和基质细胞的丰度；（2）关于免疫细胞，仅能计算一个总的免疫细胞评分，而无法给出每种免疫细胞的具体比例。
# 
# quanTIseq：(肿瘤可用)是用于根据人类 RNA-seq 数据量化肿瘤免疫状况，通过反卷积量化样本中存在的十种不同免疫细胞类型（见图1）的比例以及其他未表征细胞的比例。
# 
# TIMER：(肿瘤可用)用反卷积方法估算 32 种癌症中 6 种免疫细胞的丰度（B 细胞，CD4+ T 细胞，CD8+ T 细胞，嗜中性粒细胞，巨噬细胞和树突状细胞）。
# 
# EPIC：使用约束最小二乘回归将非负性约束条件纳入反卷积问题，从大量肿瘤基因表达数据中估算免疫和癌细胞的比例，每个样本中所有细胞分数的总和不超过一。
# 
# MCPcounter：一种基于标记基因定量肿瘤浸润免疫细胞（CD3 + T细胞，CD8 + T细胞，细胞毒性淋巴细胞，NK细胞，B淋巴细胞，单核细胞来源的细胞(单核细胞系)，髓样树突状细胞，中性粒细胞）、成纤维细胞和上皮细胞的方法。对于每个细胞类型和样本，丰度得分为细胞类型特异性基因表达值的几何平均值，独立地对每个样本进行计算的。由于分数是用任意单位表示的，它们不能直接解释为细胞分数，也不能在细胞类型之间进行比较。进行定量验证时，估计分数与真实细胞分数之间具有很高的相关性，证明了MCP-counter用于样本间比较的价值。MCP-counter已应用于对32个非血液学肿瘤的19,000多个样本中的免疫细胞和非免疫细胞进行量化。
# 
# 'CIBERSORT'、'svr' 和 'xCell' 方法对array有优化。建议array使用这三种方式

ydx_analysis_cibersort <- function(pre_data_immu=exp_data, arrays = TRUE, perm = 1000) {
  library(IOBR)
  library(tidyverse)
  library(tidyHeatmap)
  library(maftools)
  library(ggpubr)
  library(ggplot2)
  library(survival)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggsci)
  library(EPIC)
  rownames(pre_data_immu) <- toupper(rownames(pre_data_immu))
  if(arrays){
    # 将行名转换为大写字母,要不匹配不到数据
  cat("当前非肿瘤的芯片数据，正在采用芯片数据分析方法\n")
    Cibersort <- deconvo_tme(pre_data_immu, method = "cibersort",reference="lm22", arrays = TRUE, perm = perm,tumor =FALSE )
    Cibersort_abs <- deconvo_tme(pre_data_immu, method = "cibersort_abs",reference="lm22", arrays = TRUE, perm = perm,tumor =FALSE )
    xCell <- deconvo_tme(pre_data_immu, method = "xcell", arrays = TRUE,perm = perm,tumor =FALSE )
    IPS<- deconvo_tme(pre_data_immu, method = "ips", arrays = TRUE,perm = perm,tumor =FALSE )
    MCPcounter<- deconvo_tme(pre_data_immu, method = "mcpcounter", arrays = TRUE,perm = perm,tumor =FALSE ) 
    }else{
      cat("当前非肿瘤的RNA-Seq数据，正在采用RNA-Seq数据分析方法\n")
    Cibersort <- deconvo_tme(pre_data_immu, method = "cibersort",reference="lm22",arrays = FALSE,  perm = perm ,tumor =FALSE)
    Cibersort_abs<- deconvo_tme(pre_data_immu, method = "cibersort_abs",reference="lm22", arrays = FALSE,perm = perm,tumor =FALSE )
    xCell <- deconvo_tme(pre_data_immu, method = "xcell", arrays = FALSE,perm = perm ,tumor =FALSE)
    IPS<- deconvo_tme(pre_data_immu, method = "ips", arrays = FALSE,perm = perm ,tumor =FALSE)
    MCPcounter<- deconvo_tme(pre_data_immu, method = "mcpcounter", arrays = FALSE,perm = perm ,tumor =FALSE ) 
    }
  # 
  # library(RColorBrewer)
  # library(circlize)
  # library(gplots)
  # library(viridis)
  # library(oompaBase)
  # standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  #   outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  #   if (!is.null(halfwidth)) {
  #     outdata[outdata>halfwidth]=halfwidth
  #     outdata[outdata<(-halfwidth)]= -halfwidth
  #   }
  #   return(outdata)
  # }
  # 
  # 
  # cibersort=immu_analysis$Cibersort
  # MCPcounter=immu_analysis$xCell
  # xCell=immu_analysis$Cibersort_abs
  # ips=immu_analysis$IPS
  # 
  # 
  # tme_combine<-cibersort %>% 
  #   inner_join(.,MCPcounter,by       = "ID") %>% 
  #   inner_join(.,xCell,by     = "ID") %>%
  #   inner_join(.,ips,by       = "ID")
  # 
  # 

  
  return(list(Cibersort=Cibersort,Cibersort_abs=Cibersort_abs,xCell=xCell,IPS=IPS,MCPcounter=MCPcounter))
}
ydx_draw_cibersort  <- function(Cibersort=immu_analysis$Cibersort_abs,picDir = "xCell", pheatmap_width = 10,pheatmap_height =6,boxplot_width = 16,boxplot_height =8){
  library(IOBR)
  library(tidyverse)
  library(tidyHeatmap)
  library(maftools)
  library(ggpubr)
  library(ggplot2)
  library(survival)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggsci)
  library(EPIC)
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {dir.create(picDir)}
  
  
  if(length(Cibersort) == 27){
    cat("当前使用的是Cibersort_abs\n")
    Cibersort= Cibersort[,-27]
    pheatmap_name = "/CIBERSORT_abs-热图.pdf"
  }else if(length(Cibersort) == 26){
    cat("当前使用的是Cibersort\n")
    pheatmap_name = "/CIBERSORT-热图.pdf"
  }else if(length(Cibersort) == 68){
  cat("当前使用的是xCell\n")
  pheatmap_name = "/xcell-热图.pdf"  }

  
  
  cat("计算多个基因和免疫浸润的相关性\n")
  cat("此函数用来计算热图，相关性热图，分组比较以及免疫条形图\n")
  
  Cibersortname= column_to_rownames(Cibersort,"ID")
  conNum <- length(str_which(rownames(Cibersortname), "_con")) 
  cat("正常组数量:", conNum, "\n")
  treatNum <- length(str_which(rownames(Cibersortname), "_treat")) 
  cat("实验组数量:", treatNum, "\n")
  group=c(rep("con",conNum),rep("treat",treatNum))
 

  
  
  if(length(Cibersort)==68){
    re <-column_to_rownames(Cibersort[,-(66:68)],"ID")
    colnames(re) <- gsub("_xCell", "", colnames(re))
    cat("xcell图片输出\n")
  }else{
    #仅保留pvalue小于0.05的样本
    Cibersort_Pvalue <- Cibersort[Cibersort[, "P-value_CIBERSORT"] < 0.05, ]
    re <-column_to_rownames(Cibersort_Pvalue[,-(24:26)],"ID")
    colnames(re) <- gsub("_CIBERSORT", "", colnames(re))
    cat("cibersort图片输出\n")
  }

  
  
  
  re[] <- lapply(re, function(x) as.numeric(as.character(x)))
  re2 <- re[,colSums(re) > 0]   #去除丰度全为0的细胞
  
  re2 <- as.data.frame(t(re2))
  an = data.frame(group = group,
                  row.names = colnames(re2))
  
  library(pheatmap)
  
  #########################热图####################### 
  
  annotation_colors=list(Group = c("Con" = "#00C094", "Treat" = "#f58220"))
 
  pdf(paste0(picDir, pheatmap_name),width = pheatmap_width,height = pheatmap_height)
  plot1=pheatmap::pheatmap(as.matrix(re2), show_rownames = T,scale = "row",
                           annotation_col = an,
                           annotation_colors=annotation_colors,
                           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                           gaps_col = c(conNum),
                           fontsize = 10,
                           cluster_cols = F,
                           cluster_rows = T,
                           show_colnames = F,
                           cutree_rows =2)
  plot2=pheatmap::pheatmap(as.matrix(re2), show_rownames = T,scale = "row",
                           annotation_col = an,
                           annotation_colors=annotation_colors,
                           gaps_col = c(conNum),
                           fontsize = 10,
                           cluster_cols = F,
                           cluster_rows = T,
                           show_colnames = F,
                           cutree_rows =2)
  print(plot1)
  print(plot2)
  dev.off()
  
  #########################条形图绘制####################### 
  if(length(Cibersort)==68){
    re <-column_to_rownames(Cibersort[,-(66:68)],"ID")
    colnames(re) <- gsub("_xCell", "", colnames(re))
    cat("xcell图片输出\n")
  }else{
    #仅保留pvalue小于0.05的样本
    Cibersort_Pvalue <- Cibersort[Cibersort[, "P-value_CIBERSORT"] < 0.05, ]
    re <-column_to_rownames(Cibersort_Pvalue[,-(24:26)],"ID")
    colnames(re) <- gsub("_CIBERSORT", "", colnames(re))
    cat("cibersort图片输出\n")
  }
  
  
  
  
  re[] <- lapply(re, function(x) as.numeric(as.character(x)))
#  re <- re[,colSums(re) > 0]   #去除丰度全为0的细胞
  
  res <- as.data.frame(re)
  res$group = group
  res <- gather(res, key = "cell.type", value = "value",-group)
  

  
  pdf(paste0(picDir, "/CIBERSORT-分组条形图(包含未检测到的细胞类型).pdf"), height = boxplot_height, width = boxplot_width)
  
  # res <-column_to_rownames(cibersort[,-(24:26)],"ID") 
 

  plot1= ggplot(res,aes(cell.type,value,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "Cell Type", y = "Estimated Proportion") +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 10,face = "italic",colour = 'black'),
          axis.text.y = element_text(face = "italic",size = 10,colour = 'black'))+
    # 设置填充颜色的调色板
    scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
    # 添加统计分析结果，显示不同组别之间的显著差异，使用Kruskal-Wallis检验方法。
    stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")
  print(plot1)
  
  library(ggpubr)
  library(RColorBrewer)
  mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
  plot2=  ggplot(res,aes(cell.type,value,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "Cell Type", y = "Estimated Proportion") +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
    scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = group,label = ..p.signif..),method = "kruskal.test")
  print(plot2)
  
  p3 <- ggplot(res,aes(cell.type,value,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    scale_fill_manual(values = palette1[c(2,4)])+ 
    theme_bw() + 
    labs(x = "Cell Type", y = "Estimated Proportion") +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle=45,hjust = 1),
          axis.text = element_text(color = "black",size = 12))+
    stat_compare_means(aes(group = group,label = ..p.signif..),
                       method = "kruskal.test",label.y = 0.4)
  print(p3)
  
  
  
  if(length(Cibersort)== 68){
    res <-re
    res$group = group
    res <- gather(res, key = "cell.type", value = "value",-group)
    cat("xcell图片输出\n")
  }else if(length(Cibersort) == 27 | length(Cibersort) == 26){
    ##去除没有检测到的免疫细胞，转换为长数据，绘图
    re <- re[,colSums(re) > 0]   #去除丰度全为0的细胞
    res <- as.data.frame(re)
    res$group = group
    res <- gather(res, key = "cell.type", value = "value",-group)
    cat("去除未检测到的免疫细胞后画图\n")
  }
  
  plot3= ggplot(res,aes(cell.type,value,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "Cell Type", y = "Estimated Proportion") +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 10,face = "italic",colour = 'black'),
          axis.text.y = element_text(face = "italic",size = 10,colour = 'black'))+
    # 设置填充颜色的调色板
    scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
    # 添加统计分析结果，显示不同组别之间的显著差异，使用Kruskal-Wallis检验方法。
    stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")
  print(plot3)
  library(ggpubr)
  library(RColorBrewer)
  mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
  plot4=  ggplot(res,aes(cell.type,value,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "Cell Type", y = "Estimated Proportion") +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
    scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
    stat_compare_means(aes(group = group,label = ..p.signif..),method = "kruskal.test")
  print(plot4)
  
  
  library(ggpubr)
  library(stringr)
  
  p3 <- ggplot(res,aes(cell.type,value,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    scale_fill_manual(values = palette1[c(2,4)])+ 
    theme_bw() + 
    labs(x = "Cell Type", y = "Estimated Proportion") +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle=45,hjust = 1),
          axis.text = element_text(color = "black",size = 12))+
    stat_compare_means(aes(group = group,label = ..p.signif..),
                       method = "kruskal.test",label.y = 0.4)
  print(p3)
  

  dev.off()
  

  
  ####展示出每个样本的免疫细胞比例
  
  ######################################################
  ###################相关性热图#########################
  ######################################################
  ########data是保留有差异的免疫细胞的矩阵################
  # 将数据转置并转换为数据框
  if(length(Cibersort)== 68){
   
    cat("xcell暂时不输出免疫细胞百分比\n")
  }else if(length(Cibersort) == 27 | length(Cibersort) == 26){
    cat("输出免疫细胞百分比\n")
    data_corrplot <- as.data.frame(t(Cibersort_Pvalue))
    
    res <-Cibersort_Pvalue[,-(24:26)]
    name = res$ID
    colnames(res) <- gsub("_CIBERSORT", "", colnames(res))
    res <-res[,-1]
    res[] <- lapply(res, function(x) as.numeric(as.character(x)))
    ciber.res <- res[,colSums(res) > 0]   #去除丰度全为0的细胞
    library(corrplot)
    
    
    #install.packages('corrplot')
    library(corrplot)
    # 创建 PDF 设备并设置大小
    pdf(paste0(picDir, "/CIBERSORT_corHeatmap_all_所有相关性热图.pdf"), height = 8, width = 8)
    # 定义methods向量
    methods <- c("circle", "pie", "square", "ellipse", "number", "shade", "color")
    # 打开PDF设备，指定文件名
    for (method in methods) {
      # 绘制相关性热图
      plot1 = corrplot(corr=round(cor(ciber.res),2),
                       method = method,
                       order = "hclust",
                       tl.col = "black",
                       addCoef.col = "black",
                       number.cex = 0.5,
                       tl.cex = 0.5,
                       col = colorRampPalette(c("navy", "white", "firebrick3"))(50))
      
      plot2 = corrplot.mixed(corr=round(cor(ciber.res),2),# 计算相关性矩阵并保留两位小数
                             lower.col = "black", #左下方字体颜色为黑色
                             tl.pos = "lt",  #标签出现在左侧和顶部
                             number.cex = 0.5, #左下方字号为0.5
                             upper = method,
                             tl.cex = 0.5) #标签字号为0.5
      print(plot1)
      print(plot2)
    }
    # 结束绘图，关闭PDF设备
    dev.off()
    cat("结束绘图，关闭PDF设备\n")
    
    ######################################比例图
    Cibersort_Pvalue <- Cibersort[Cibersort[, "P-value_CIBERSORT"] < 0.05, ]
    # 首先变为长数据
    cibersort_long <- Cibersort_Pvalue %>% 
      select(`P-value_CIBERSORT`,Correlation_CIBERSORT, RMSE_CIBERSORT,ID,everything()) %>% 
      pivot_longer(- c(1:4),names_to = "cell_type",values_to = "fraction") %>% 
      dplyr::mutate(cell_type = gsub("_CIBERSORT","",cell_type),
                    cell_type = gsub("_"," ",cell_type))
    
    pdf(paste0(picDir, "/CIBERSORT_免疫细胞比例图.pdf"), height = 8, width = 10)
    p1 <- cibersort_long %>% 
      ggplot(aes(ID,fraction))+
      geom_bar(stat = "identity",position = "stack",aes(fill=cell_type))+
      labs(x=NULL)+
      scale_y_continuous(expand = c(0,0))+
      scale_fill_manual(values = palette4,name=NULL)+ # iobr还给大家准备了几个色盘，贴心！
      theme_bw()+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "bottom"
      )
    print(p1)
    p1 <- cell_bar_plot(input = Cibersort_Pvalue,  palette = 1,legend.position = "bottom",title = "CIBERSORT Cell Fraction")
    p2 <- cell_bar_plot(input = Cibersort_Pvalue,  palette = 2,legend.position = "bottom",title = "CIBERSORT Cell Fraction")
    p3 <- cell_bar_plot(input = Cibersort_Pvalue,  palette = 3,legend.position = "bottom",title = "CIBERSORT Cell Fraction")
    p4 <- cell_bar_plot(input = Cibersort_Pvalue,  palette = 4,legend.position = "bottom",title = "CIBERSORT Cell Fraction")
    
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    dev.off()
  }
  
 
  
}


ydx_draw_cibersort_cor <- function(pre_data_immu=exp_data,
                                   Cibersort= immu_analysis$Cibersort,
                                   genename=t1dm_mci_wgcna_mit,
                                   picDir = "免疫浸润相关性分析",
                                   height = 8, width = 16,box_width=10,box_height=7){
  
  if(length(Cibersort) == 27){
    cat("当前使用的是Cibersort_abs")
    Cibersort= Cibersort[,-27]
  }
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {dir.create(picDir)}
  
  Cibersortname= column_to_rownames(Cibersort,"ID")
  conNum <- length(str_which(rownames(Cibersortname), "_con")) 
  cat("正常组数量:", conNum, "\n")
  treatNum <- length(str_which(rownames(Cibersortname), "_treat")) 
  cat("实验组数量:", treatNum, "\n")
  group=c(rep("con",conNum),rep("treat",treatNum))

  #仅保留pvalue小于0.05的样本
  Cibersort_Pvalue <- Cibersort[Cibersort[, "P-value_CIBERSORT"] < 0.05, ]
  
  tcga_gsva <-column_to_rownames(Cibersort_Pvalue[,-(24:26)],"ID")
  tcga_gsva[] <- lapply(tcga_gsva, function(x) as.numeric(as.character(x)))
  colnames(tcga_gsva) <- gsub("_CIBERSORT", "", colnames(tcga_gsva))
 
  
  cat("仅保留实验组的数据做免疫细胞相关性分析\n")
  treat_rows <- tcga_gsva[grep("treat", rownames(tcga_gsva)), ]
  pre_data_immu_treat  <- pre_data_immu[,grep("treat", colnames(pre_data_immu)) ]
  tcga_gsva <- treat_rows[,colSums(treat_rows) > 0]   #去除丰度全为0的细胞
  
  
  
  
  pre_data_immu_treat <- as.data.frame(t(pre_data_immu_treat[rownames(pre_data_immu_treat) %in% genename,]))
  
  


  identical(rownames(tcga_gsva),rownames(pre_data_immu_treat))
  

  
  cor_res <- linkET::correlate(pre_data_immu_treat, tcga_gsva,method = "spearman")
  
p =   linkET::qcorrplot(cor_res) +
  linkET::geom_square() +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))
  

df_r <- cor_res$r %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1,names_to = "cell_type",values_to = "correlation")
df_p <- cor_res$p %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1,names_to = "cell_type",values_to = "pvalue")

df_cor <- df_r %>% 
  left_join(df_p) %>% 
  mutate(stars = cut(pvalue,breaks = c(-Inf,0.05,0.01,0.001,Inf),right = F,labels = c("***","**","*"," ")))
## Joining with `by = join_by(gene, cell_type)`
library(ggplot2)
pdf(paste0(picDir, "/疾病组hub基因与免疫细胞的相关性.pdf"), height = height, width = width)
p1 =  ggplot(df_cor, aes(cell_type,gene))+
  geom_tile(aes(fill=correlation))+
  geom_text(aes(label=stars), color="black", size=4)+
  scale_fill_gradient2(low='#67B26F', high='#F2AA9D',mid = 'white',
                       limit=c(-1,1),name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_text(size=8,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())
  

  library(ggplot2)
  library(dplyr)
  p2 = ggplot(df_cor, aes(cell_type, gene)) + 
    geom_tile(aes(fill = correlation), colour = "white",size=1)+
    scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
    geom_text(aes(label=stars),col ="black",size = 5)+
    theme_minimal()+# 不要背景
    theme(axis.title.x=element_blank(),#不要title
          axis.ticks.x=element_blank(),#不要x轴
          axis.title.y=element_blank(),#不要y轴
          axis.text.x = element_text(angle = 45, hjust = 1),# 调整x轴文字
          axis.text.y = element_text(size = 8))+#调整y轴文字
    #调整legen
    labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))

# 绘制热图
p3 <- ggplot(df_cor, aes(cell_type, gene)) + 
  geom_tile(aes(fill = correlation), colour = "black", size = 1) +
  scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e41a1c") +
  geom_text(aes(label = stars), col = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(fill = "Correlation", subtitle = paste0(" * p < 0.05","\n\n","** p < 0.01"))
  
  library(ggplot2)
  library(dplyr)
  p4 = ggplot(df_cor, aes(cell_type, gene)) + 
    geom_tile(aes(fill = correlation),size=1)+
    scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
    geom_text(aes(label=stars),col ="black",size = 5)+
    theme_minimal()+# 不要背景
    theme(axis.title.x=element_blank(),#不要title
          axis.ticks.x=element_blank(),#不要x轴
          axis.title.y=element_blank(),#不要y轴
          axis.text.x = element_text(angle = 45, hjust = 1),# 调整x轴文字
          axis.text.y = element_text(size = 8))+#调整y轴文字
    #调整legend
    labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
  
  
  #################每个基因单独跑一遍
  cat("开始绘制单个基因的相关性分析 \n")
  #仅保留有统计学意义的分组pvalue小于0.05的样本
  Cibersort_Pvalue <- Cibersort[Cibersort[, "P-value_CIBERSORT"] < 0.05, ]
  tcga_gsva <-column_to_rownames(Cibersort_Pvalue[,-(24:26)],"ID")
  tcga_gsva[] <- lapply(tcga_gsva, function(x) as.numeric(as.character(x)))
  tcga_gsva <- tcga_gsva[,colSums(tcga_gsva) > 0]   #去除丰度全为0的细胞
  
  
  colnames(tcga_gsva) <- gsub("_CIBERSORT", "", colnames(tcga_gsva))
  
  

  
  
  for (gene in genename) {
    if(all(genename == "")){print("gene_name为空，后续基因调控免疫细胞的分析取消")
    }else{
      
      cat("将疾病组",gene,"分为高低组，观察免疫细胞的变化\n")
      library(limma)
      library(reshape2)
      library(ggpubr)
      library(vioplot)
      library(ggExtra)
      gene_exp <- pre_data_immu[gene, ]
      gene_exp= as.data.frame(gene_exp)
      gene_exp$type <- ifelse(grepl("_con", rownames(gene_exp)), "Normal", "Treat")
      #获取疾病组
      tumorData=gene_exp[gene_exp$type=="Treat",1,drop=F]
      
      tumorData=as.matrix(tumorData)
      data=avereps(tumorData)
      #按照预定，将基因分为高低表达组
      data=as.data.frame(data)
      name =colnames(gene_exp)[1]
      data$gene=ifelse(data[,name]>median(data[,name]), "High", "Low")
      
      treat_immune <- tcga_gsva[grep("_treat", rownames(tcga_gsva)), ]
      treat_immune =as.data.frame(treat_immune)
      treat_immune=avereps(treat_immune)
      treat_immune =as.data.frame(treat_immune)
      
      
      
      
      # 将基因表达与免疫细胞表达结合起来
      sameSample=intersect(row.names(treat_immune), row.names(data))
      rt=cbind(tcga_gsva[sameSample,,drop=F], data[sameSample,,drop=F])
      
      data=rt
      data=melt(data,id.vars=c("gene"))
      colnames(data)=c("gene", "Immune", "Expression")

      
      
      
      
      picDir_dir = "单个基因免疫浸润分析"
      # 拼接子目录路径
      picDir_dir_full <- file.path(picDir, picDir_dir)
      # 检查 picDir_dir 子目录是否存在，如果不存在则创建
      if (!file.exists(picDir_dir_full)) {
        dir.create(picDir_dir_full)
      }

      
      outTab=data.frame()
      for(i in colnames(rt)[1:(ncol(rt)-2)]){
        x=as.numeric(rt[,name])
        y=as.numeric(rt[,i])
        if(sd(y)==0){y[1]=0.00001}
        cor=cor.test(x, y, method="spearman")
        outVector=cbind(Cell=i, cor=cor$estimate, pvalue=cor$p.value)
        outTab=rbind(outTab,outVector)
        if(cor$p.value<0.05){
          outFile=paste0(picDir_dir_full,"/",gene,"_cor.", i, "相关性",round(cor$estimate,2),".pdf")
          df1=as.data.frame(cbind(x,y))
          p1=ggplot(df1, aes(x, y)) + 
            xlab(paste0(gene, " expression")) + ylab(i)+
            geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
            stat_cor(method = 'spearman', aes(x =x, y =y))
          p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
          pdf(file=outFile, width=5.2, height=5)
          print(p2)
          dev.off()
        }
      }
      outTab
      rownames(outTab) <- NULL
   
 
      
  
      
      contains_gene_exp <- grepl("gene_exp", data$Immune)
      
      data3=data
      data2 <- data3[!contains_gene_exp, ]
      
      
      library(ggpubr)
      #绘制boxplot
      pdf(paste0(picDir,"/",gene,"----基因高低分组比较.pdf"), width=box_width, height=box_height)  
      p=ggboxplot(data2, x="Immune", y="Expression", fill = "gene", 
                  ylab="Infiltration levels",
                  xlab="",
                  legend.title=gene,
                  palette = c("#0072B5CC","#E18727CC"),
                  add = "none")
      p=p+rotate_x_text(60)
      p1 <- p + ggtitle("Expression of Immune Genes") +
        ggpubr::stat_compare_means(aes(group = gene),
                                   method = "wilcox.test",
                                   symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                                   label = "p.signif")
      print(p1)

      
      
      
      plot1= ggplot(data2,aes(Immune,Expression,fill = gene)) + 
        geom_boxplot(outlier.shape = 21,color = "black") + 
        theme_bw() + 
        labs(x = "Cell Type", y = "Estimated Proportion") +
        theme(legend.position = "top") + 
        theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 10,face = "italic",colour = 'black'),
              axis.text.y = element_text(face = "italic",size = 10,colour = 'black'))+
        # 设置填充颜色的调色板
        scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
        # 添加统计分析结果，显示不同组别之间的显著差异，使用Kruskal-Wallis检验方法。
        stat_compare_means(aes(group = gene),label = "p.format",size=3,method = "wilcox.test")
      
      print(plot1)
      
      library(ggpubr)
      library(RColorBrewer)
      mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
      plot2=  ggplot(data2,aes(Immune,Expression,fill = gene)) + 
        geom_boxplot(outlier.shape = 21,color = "black") + 
        theme_bw() + 
        labs(x = "Cell Type", y = "Estimated Proportion") +
        theme(legend.position = "top") + 
        theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
        scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = gene,label = ..p.signif..),method = "kruskal.test")
      print(plot2)
      
      p3 <- ggplot(data2,aes(Immune,Expression,fill = gene)) + 
        geom_boxplot(outlier.shape = 21,color = "black") + 
        scale_fill_manual(values = palette1[c(2,4)])+ 
        theme_bw() + 
        labs(x = "Cell Type", y = "Estimated Proportion") +
        theme(legend.position = "top") + 
        theme(axis.text.x = element_text(angle=45,hjust = 1),
              axis.text = element_text(color = "black",size = 12))+
        stat_compare_means(aes(group = gene,label = ..p.signif..),
                           method = "kruskal.test",label.y = 0.4)
      print(p3)
      
      
      
      
      
      
      dev.off()
      
      
 
      
      
      
      
      
      
      
      
      
      

      #读取输入文件
      data = outTab
      data$cor <- as.numeric(data$cor)
      data$pvalue <- as.numeric(data$pvalue)
      #定义圆圈颜色的函数
      p.col = RColorBrewer::brewer.pal(5,"Set3")
      fcolor = function(x,p.col){
        color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                                   ifelse(x>0.2,p.col[4], p.col[5])
        )))
        return(color)
      }
      
      #定义设置圆圈大小的函数
      p.cex = seq(2.5, 5.5, length=5)
      fcex = function(x){
        x=abs(x)
        cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                                 ifelse(x<0.4,p.cex[4],p.cex[5]))))
        return(cex)
      }
      
      #根据pvalue定义圆圈的颜色
      points.color = fcolor(x=data$pvalue,p.col=p.col)
      data$points.color = points.color
      
      #根据相关系数定义圆圈的大小
      points.cex = fcex(x=data$cor)
      data$points.cex = points.cex
      data=data[order(data$cor),]
      
      ########绘制图形########
      xlim = ceiling(max(abs(data$cor))*10)/10         #x轴范围    
      #输出图形
      
      pdf(paste0(picDir,"/",gene,"----棒棒糖图.pdf"), width=9, height=7)  
      layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
      par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
      plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
      rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
      grid(ny=nrow(data),col="white",lty=1,lwd=2)
      #绘制图形的线段
      segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
      #绘制图形的圆圈
      points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
      #展示免疫细胞的名称
      text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5)
      #展示pvalue
      pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
      redcutoff_cor=0
      redcutoff_pvalue=0.05
      text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
      axis(1,tick=F)
      
      #绘制圆圈大小的图例
      par(mar=c(0,4,3,4))
      plot(1,type="n",axes=F,xlab="",ylab="")
      legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")
      
      #绘制圆圈颜色的图例
      par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
      barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
      axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
      dev.off()
      

      
      
      
    }
    
  }
  }
  
  
  
  
  
  
  












  







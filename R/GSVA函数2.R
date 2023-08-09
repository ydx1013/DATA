ydx_GSVA  <-  function(exp_data,picDir = "GSVA",category  = "kegg",qianzui= "KEGG_",width = 6, height = 8,adjpvalue_cut  = 0.05,logFC_cut = 0){
  library(msigdbr)
  library(dplyr)
  library(data.table)
  library(GSVA)
  library(limma)
  library(stringr)
  library(ggplot2)
  library(easyTCGA)
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }

  category <- toupper(category)
  switch(category,
         "KEGG" = {
           category = "C2"
           subcategory = "CP:KEGG"
           
         },
         "BP" = {
           category = "C5"
           subcategory = "GO:BP"
           qianzui = "GOBP_"
         },
         "CC" = {
           category = "C5"
           subcategory = "GO:CC"
           qianzui = "GOCC_"
         },
         "MF" = {
           category = "C5"
           subcategory = "GO:MF"
           qianzui = "GOMF_"
         },
         "immu" = {
           category = "C7"
           subcategory = "IMMUNESIGDB"
           qianzui = ""
         },
         "REACTOME" = {
           category = "C2"
           subcategory = "CP:REACTOME"
           qianzui = "REACTOME_"
         },
         "C5" = {
           category = "C5"
           subcategory = "NULL"
           qianzui = ""
         }

  )
  

  GeneName = as.character(rownames(exp_data))
  result <- sum(grepl("[a-z]", GeneName))
  if(result>length(GeneName)/2){
    species = "Mus musculus"
    cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
  }else{
    species = "Homo sapiens"
    cat("检测到基因属于人类，使用人类基因集做分析\n")
  }
  exp_data <- as.matrix(exp_data)
  exp_data <- avereps(exp_data)
  exp_data <- normalizeBetweenArrays(exp_data)

  #查看msigdbr包里自带的物种
  msigdbr_species()
  #其他物种可自行到GSEA官网下载，http://software.broadinstitute.org/gsea/downloads.jsp
  #######生成KEGG  GENESET#######################
  subcategory = 
    if (subcategory == "NULL") {
      h <- msigdbr(species = species,# 物种拉丁名
                   category = category)
    } else {
      h <- msigdbr(species = species, category = category, subcategory = subcategory)
    }
  
  
  # 示例数据表达矩阵的基因名是gene symbol，这里就选gene_symbol。
  # 如果你的表达矩阵以ENTREZ ID作为基因名，就把下面这段的gene_symbol换成entrez_gene
  h <- select(h, gs_name, gene_symbol) %>% #或entrez_gene
    as.data.frame %>% 
    split(., .$gs_name) %>% 
    lapply(., function(x)(x$gene_symbol)) #或entrez_gene
  
  # 在每个geneset里面去掉重复的基因
  gs <- lapply(h, unique)
  
  # 接下来去掉那些在两个或更多个pathways里出现过的genes
  count <- table(unlist(gs))
  keep <- names(which(table(unlist(gs)) < 2))
  gs <- lapply(gs, function(x) intersect(keep, x))
  
  # 过滤之后，很多pathway一个gene都不剩了，去掉这些
  gs <- gs[lapply(gs, length) > 0]
  
  # 这一句就完成了GSVA分析
  gsva_es <- gsva(as.matrix(exp_data), gs)
  # 分组
  
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
  group=c(rep("con",conNum),rep("treat",treatNum))
  Type=as.factor(group)
  
  DEG_list <- diff_analysis(exprset = gsva_es,
                            group = Type,
                            is_count = FALSE,  # 不是 count 数据
                            logFC_cut = logFC_cut,  # 可以直接筛选结果
                            adjpvalue_cut = adjpvalue_cut)
  
  # 假设DEG_list是包含多个数据框的列表
  # 使用gsub()函数去除DEG_list中第一个数据框中所有的"kegg_"字符串
  DEG_list[[1]][] <- lapply(DEG_list[[1]], function(x) gsub(qianzui, "", x))
  DEG_list[[2]][] <- lapply(DEG_list[[2]], function(x) gsub(qianzui, "", x))
  
  #把通路的limma分析结果保存到文件
  write.csv(gsva_es, paste0(picDir,"/gsva_es打分矩阵.csv"), quote = F)
  write.csv(DEG_list$deg_limma, paste0(picDir,"/gsva_limma.csv"), quote = F)
  write.csv(DEG_list$deg_wilcoxon, paste0(picDir,"/gsva_wilcoxon.csv"), quote = F)
  
  
  #### 对GSVA的差异分析结果进行热图可视化 #### 
  
  keep <-rownames(DEG_list$deg_limma)
  
  dat <- gsva_es[keep,] #选取前50进行展示

  rownames(dat) <- gsub(qianzui, "", rownames(dat))
  
  library("pheatmap")
  annotation_col= data.frame(Group = factor(c(rep("Con",conNum),rep("Treat",treatNum))))
  row.names(annotation_col) = colnames(dat)
  annotation_colors=list(Group = c("Con" = "#00C094", "Treat" = "#f58220"))
  pheatmap::pheatmap(dat, show_rownames = T,scale = "row",
           annotation_col = annotation_col,
           annotation_colors=annotation_colors,
           gaps_col = c(conNum),
           fontsize = 10,
           cluster_cols = F,
           cluster_rows = T,
           show_colnames = F,
           cutree_rows =2,
           filename = paste0(picDir,"/pheatmap_by_GSVA_sd10.pdf"),width = 10,height = 6)
  
  

  

  
  #### 发散条形图绘制 ####
  library(tidyverse)  # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
  library(ggthemes)
  library(ggprism)
  degs <- DEG_list$deg_limma  #载入gsva的差异分析结果
  Diff <- rbind(subset(degs,logFC>0)[1:20,], subset(degs,logFC<0)[1:20,]) #选择上下调前20通路     
  dat_plot <- data.frame(id  = row.names(Diff),
                         p   = as.numeric(Diff$P.Value),# 将P列转换为数字
                         lgfc= Diff$logFC)
  
  dat_plot$group <- ifelse(dat_plot$lgfc>0 ,1,-1)    # 将上调设为组1，下调设为组-1

  dat_plot$lg_p <- -log10(dat_plot$p)*dat_plot$group # 将上调-log10p设置为正，下调-log10p设置为负
  
  dat_plot$id <- str_replace(dat_plot$id, qianzui,"");dat_plot$id[1:10]
  dat_plot$threshold <- factor(ifelse(abs(dat_plot$p) <= adjpvalue_cut,
                                      ifelse(dat_plot$lgfc >0 ,'Up','Down'),'Not'),
                               levels=c('Up','Down','Not'))
  
  dat_plot <- dat_plot %>% arrange(lg_p)
  dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
  
  ## 设置不同标签数量
  low1 <- dat_plot %>% filter(lg_p < log10(adjpvalue_cut)) %>% nrow()
  low0 <- dat_plot %>% filter(lg_p < 0) %>% nrow()
  high0 <- dat_plot %>% filter(lg_p < -log10(adjpvalue_cut)) %>% nrow()
  high1 <- nrow(dat_plot)
  p <- ggplot(data = dat_plot,aes(x = id, y = lg_p, 
                                  fill = threshold)) +
    geom_col()+
    coord_flip() + 
    scale_fill_manual(values = c('Up'= '#36638a','Not'='#cccccc','Down'='#7bcd7b')) +
    geom_hline(yintercept = c(-log10(adjpvalue_cut),log10(adjpvalue_cut)),color = 'white',size = 0.5,lty='dashed') +
    xlab('') + 
    ylab('-log10(P.Value) of GSVA score') + 
    guides(fill="none")+
    ggprism::theme_prism(border = T) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
              hjust = 0,color = 'black') + #黑色标签
    geom_text(data = dat_plot[(low1 ):low0,],aes(x = id,y = 0.1,label = id),
              hjust = 0,color = 'black') + # 灰色标签
    geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
              hjust = 1,color = 'black') # 黑色标签

  ggsave(paste0(picDir,"/GSVA_barplot_pvalue.pdf"),p,width = 12,height  = 8)
  
return(DEG_list=DEG_list)
  }
  
ydx_friends <- function(GeneName=t1dm_mci_wgcna_mit, picDir = "GO",width= 8,height= 10,fromType = "SYMBOL",toType = "ENTREZID"){
  #加载包
  
  cat("背景资料：https://mp.weixin.qq.com/s/_Wt_GmC8yjcvEdXBNRUQTw\n")
  library(GOSemSim)
  library(reshape2)
  library(ggplot2)
  library(clusterProfiler)
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
  if(is.data.frame(GeneName)){
    GeneName=GeneName[,1]
  }
  GeneName = as.character(GeneName)
  # 判断是否包含小写字母并返回结果为TRUE的行
  result <- sum(grepl("[a-z]", GeneName))
  if(result>length(GeneName)/2){
    library(org.Mm.eg.db)
    org_db= "org.Mm.eg.db"
    cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
  }else{
    library(org.Hs.eg.db)
    org_db= "org.Hs.eg.db"
    cat("检测到基因属于人类，使用人类基因集做分析\n")
  }
  ENSEMBLlist <- bitr(GeneName, fromType = fromType, toType = toType, org_db)
  ENSEMBLlist_ENTREZID <- ENSEMBLlist$ENTREZID
  #用godata()函数来构建相应物种的Molecular Function本体的GO DATA
  mf <- godata(org_db, ont="MF", computeIC = FALSE)
  #用godata()函数来构建相应物种的Cellular Component本体的GO DATA
  cc <- godata(org_db, ont="CC", computeIC = FALSE)
  #用godata()函数来构建相应物种的Biological Process本体的GO DATA
  bp <- godata(org_db, ont="BP", computeIC = FALSE)
 
   ########计算语义相似度
  #用mgeneSim来计算MF本体，基因之间的语义相似度，结果为一个行列相同的矩阵
  simmf <- mgeneSim(ENSEMBLlist_ENTREZID, semData = mf, measure = "Wang", drop = NULL, combine = "BMA")
  #用mgeneSim来计算CC本体，基因之间的语义相似度，结果为一个行列相同的矩阵
  simcc <- mgeneSim(ENSEMBLlist_ENTREZID, semData = cc, measure = "Wang", drop = NULL, combine = "BMA")
  #用mgeneSim来计算BP本体，基因之间的语义相似度，结果为一个行列相同的矩阵
  simbp <- mgeneSim(ENSEMBLlist_ENTREZID, semData = bp, measure = "Wang", drop = NULL, combine = "BMA")
  
  cat("选择共有的基因做分析\n")
  # 使用Reduce()和intersect()获取共有的列名
  common_columns <- Reduce(intersect, list(colnames(simbp), colnames(simcc), colnames(simmf)))
  
  # 筛选出共有列名对应的数据框
  simbp_filtered <- simbp[common_columns, common_columns]
  simcc_filtered <- simcc[common_columns, common_columns]
  simmf_filtered <- simmf[common_columns, common_columns]
  
  
  #或者计算基因在MF、CC、BP本体下的几何平均值
  fsim <- (simmf_filtered * simcc_filtered * simbp_filtered)^(1/3)
  
  # 1. 将fsim的行名和列名与ENSEMBLlist的第一列匹配，获取匹配成功的索引
  matched_rows <- match(rownames(fsim), ENSEMBLlist$ENTREZID)
  matched_cols <- match(colnames(fsim), ENSEMBLlist$ENTREZID)
  
  # 2. 使用匹配的索引将fsim中的行名和列名替换为ENSEMBLlist的第二列的数据
  rownames(fsim) <- ENSEMBLlist$SYMBOL[matched_rows]
  colnames(fsim) <- ENSEMBLlist$SYMBOL[matched_cols]
  
  fsim_data=as.data.frame(fsim)
  #将基因自己和自己的相似度设为NA，方便接下来去掉。
  for (i in 1:ncol(fsim)){
    fsim[i,i] <- NA
  }
  
  y <- melt(fsim) #把宽格式数据转化成长格式，其实就是把正方形矩阵转成三列
  y <- y[!is.na(y$value),] #删掉带NA的行
  
  
  ## 开始画图
  #计算每个基因跟其他基因相似度的平均值
  y.mean <- aggregate(.~Var1,y,mean) 
  m <- y.mean$value
  names(m) <- y.mean$Var1
  #按平均值给基因名排序，便于画图
  y$Var1 <- factor(y$Var1, levels=names(sort(m)))
  
  f <- function(y) {
    r <- quantile(y, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    r[3] <- mean(y)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  colo <- c("#FF8F00",
            "#FFA700", "#FFBF00", "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00",
            "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", "#20FF00",
            "#08FF00", "#00FF10", "#00FF28", "#00FF40", "#00FF58", "#00FF70", "#00FF87",
            "#00FF9F", "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
            "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", "#0040FF", "#0028FF",
            "#0010FF", "#0800FF", "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF",
            "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", "#FF00D7",
            "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", "#FF0060", "#FF0048", "#FF0030",
            "#FF0018")
  
  color=colo[1:23]
  
  p1 <- ggplot(y, aes(Var1, value, fill = factor(Var1))) + 
  #  ggsci::scale_fill_ucscgb() +
    guides(fill=FALSE) + #不显示图例
    stat_summary(fun.data= f, geom='boxplot') + 
    geom_hline(aes(yintercept=0.75), linetype="dashed") + #画一条虚线
    coord_flip() + # x、y坐标轴互换
    xlab("") + ylab("") + 
    theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
    theme_bw() + 
    theme(panel.border=element_rect(size=1)) #边框粗细 
  p1
  
  # 保存到pdf文件
  ggsave(paste0(picDir,"/friends_box.pdf"),width= width,height=height )
  return(list(fsim=fsim_data,table = y))
  
}



ydx_draw_GO <- function(GeneName,  qvalueFilter=1, pvalueFilter=1, picDir = "GO",fromType = "SYMBOL",toType = "ENTREZID") {
  #相关R包载入：
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
if(is.data.frame(GeneName)){
  GeneName=GeneName[,1]
}
    
  #GeneName = mit_gene
  GeneName = as.character(GeneName)
  packages <- c("topGO", "ReactomePA", "tidyverse", "clusterProfiler", "biomaRt", "enrichplot", "dplyr", "data.table", "ggplot2", "patchwork", "cols4all", "ggraph", "igraph","RColorBrewer", "ComplexHeatmap", "ggpubr", "circlize", "org.Hs.eg.db", "org.Mm.eg.db")
  lapply(packages, require, character.only = TRUE)
  # 转换基因名称格式
  # 判断是否包含小写字母并返回结果为TRUE的行
  result <- sum(grepl("[a-z]", GeneName))
  if(result>length(GeneName)/2){
    org_db= "org.Mm.eg.db"
    cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
  }else{
    org_db= "org.Hs.eg.db"
    cat("检测到基因属于人类，使用人类基因集做分析\n")
  }
  ENSEMBLlist <- bitr(GeneName, fromType = fromType, toType = toType, org_db)
  ENSEMBLlist <- ENSEMBLlist$ENTREZID
  

  # 进行GO分析
  ego <- enrichGO(ENSEMBLlist,ont="ALL",OrgDb = org_db)
  ego_MF <-enrichGO(ENSEMBLlist,ont="MF",OrgDb = org_db, qvalueCutoff = qvalueFilter, pvalueCutoff = pvalueFilter)
  ego_BP <-enrichGO(ENSEMBLlist,ont="BP",OrgDb = org_db, qvalueCutoff = qvalueFilter, pvalueCutoff = pvalueFilter)
  ego_CC <-enrichGO(ENSEMBLlist,ont="CC",OrgDb = org_db, qvalueCutoff = qvalueFilter, pvalueCutoff = pvalueFilter)


  
  # 进行KEGG分析
  if (org_db == "org.Hs.eg.db") {
    kegg <- enrichKEGG(ENSEMBLlist, organism = "hsa", qvalueCutoff=qvalueFilter,pvalueCutoff=pvalueFilter, keyType = "kegg", pAdjustMethod = "BH")
  } else if (org_db == "org.Mm.eg.db") {
    kegg <- enrichKEGG(ENSEMBLlist, organism = "mmu", qvalueCutoff=qvalueFilter,pvalueCutoff=pvalueFilter, keyType = "kegg", pAdjustMethod = "BH")
    # 使用dplyr的mutate()和str_replace()函数去除指定值
    kegg@result$Description <- str_replace_all(kegg@result$Description, "- Mus musculus \\([^\\)]+\\)", "")    
    } else {
    stop("Not a valid database")
  }
  

 
  
  # 将富集分析结果写入Excel文件
  ego_data_frames <- list("ego" = ego, "ego_MF" = ego_MF, "ego_CC" = ego_CC, "ego_BP" = ego_BP,"kegg"=kegg)
  for (df_name in names(ego_data_frames)) {
    df <- as.data.frame(ego_data_frames[[df_name]])
    write.table(df, paste0(picDir, "/", df_name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  ego@result <- ego@result[order(ego@result$p.adj), ]
  
  # 生成散点图和条形图以可视化富集分析结果
  ego_list <- list(ego, ego_MF, ego_BP, ego_CC, kegg)
  title_list <- c("EnrichmentGO_ALL_dot", "EnrichmentGO_MF_dot", "EnrichmentGO_BP_dot", "EnrichmentGO_CC_dot","EnrichmentGO_kegg_dot")
  file_list <- paste0(picDir, "/", title_list, ".pdf")   # 修改保存文件路径
  width <- 8
  height <- 7
  
  for(i in 1:length(ego_list)){
   # p <- dotplot(ego_list[[i]], title=title_list[i], orderBy = "x", label_format = 100)
   # ggsave(p, file=file_list[i], width=width, height=height)
    p <- barplot(ego_list[[i]], showCategory=10, label_format = 100, title=title_list[i])
    ggsave(p, file=paste0(file_list[i], "_bar.pdf"), width=width, height=height)
  }

  
  
  #######################################################################################
  #######################################################################################
  ############################# 绘制GO圈图#############################
  ####################################################################################### 
  #######################################################################################
  ontology.col=c("#00AFBB", "#E7B800", "#90EE90")
  data=ego[order(ego$p.adjust),]
  datasig=data[data$p.adjust<0.05,,drop=F]
  BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
  CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
  MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
  BP = head(BP,6)
  CC = head(CC,6)
  MF = head(MF,6)
  data = rbind(BP,CC,MF)
  main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]
  
  #整理圈图数据
  BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
  Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
  ratio = Gene/BgGene
  logpvalue = -log(data$pvalue,10)
  logpvalue.col = brewer.pal(n = 8, name = "Reds")
  f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
  BgGene.col = f(logpvalue)
  df = data.frame(ego=data$ID,start=1,end=max(BgGene))
  rownames(df) = df$ego
  bed2 = data.frame(ego=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
  bed3 = data.frame(ego=data$ID,start=1,end=Gene,BgGene=Gene)
  bed4 = data.frame(ego=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
  bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5
  
  #绘制圈图的主体部分
  pdf(paste0(picDir,"/GO.circlize.pdf"),width=10,height=10)
  par(omi=c(0.1,0.1,0.1,1.5))
  circos.par(track.margin=c(0.01,0.01))
  circos.genomicInitialize(df,plotType="none")
  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
  }, track.height = 0.08, bg.border = NA,bg.col = main.col)
  
  for(si in get.all.sector.index()) {
    circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
                major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
  }
  f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
  circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                           border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                      })
  circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                           border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                      })
  circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                      panel.fun = function(region, value, ...) {
                        cell.xlim = get.cell.meta.data("cell.xlim")
                        cell.ylim = get.cell.meta.data("cell.ylim")
                        for(j in 1:9) {
                          y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                          circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                        }
                        circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                           border = NA, ...)
                        #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                      })
  circos.clear()
  #绘制圈图中间的图例
  middle.legend = Legend(
    labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
    type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
    title="",nrow=3,size= unit(3, "mm")
  )
  circle_size = unit(1, "snpc")
  draw(middle.legend,x=circle_size*0.42)
  #绘制GO分类的图例
  main.legend = Legend(
    labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
    legend_gp = gpar(col=ontology.col), title_position = "topcenter",
    title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
    grid_width = unit(5, "mm")
  )
  #绘制显著性pvalue的图例
  logp.legend = Legend(
    labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
    type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
    title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
    size = unit(3, "mm")
  )
  lgd = packLegend(main.legend,logp.legend)
  circle_size = unit(1, "snpc")
  print(circle_size)
  draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
  dev.off()
  #######################################################################################
  #######################################################################################
  ############################# 第一种画法#############################
  ####################################################################################### 
  #######################################################################################
  ego_list <- list(ego, ego_MF, ego_BP, ego_CC, kegg)
  title_list = c("pink_GO_BP enrichment barplot", "pink_GO_CC enrichment barplot", "pink_ego enrichment barplot", "pink_GO_MF enrichment barplot", "pink_KEGG enrichment barplot")
for (i in 1:length(ego_list)) {
  dt <- as.data.frame(ego_list[[i]])
  # 按照 -log10(pvalue) 排序
  dt <- dt[order(dt$pvalue, decreasing = FALSE), ]
  
  # 取前 15 条通路
  dt <- head(dt, n = 15)
  
  
  #按照通路指定绘图顺序(转换为因子)：
  dt$Description <- factor(dt$Description,
                           levels = rev(dt$Description))
  #先自定义主题：
  mytheme <- theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14,
                              hjust = 0.5,
                              face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
  )
  #常规富集条形图绘图：
  p <- ggplot(data = dt, aes(x = Count, y = Description, fill = -log10(pvalue))) +
    scale_fill_distiller(palette = "RdPu",direction = 1) + #更改配色
    geom_bar(stat = "identity", width = 0.8) + #绘制条形图
    labs(x = "Number of Gene", y = "", title = title_list[i]) + #修改/添加各标题
    theme_bw() + mytheme #主题更改
  p
  
  # 保存为PDF文件
  ggsave(file = paste0(picDir, "/",title_list[i], ".pdf"), plot = p, width = 10, height = 6, units = "in")
}
  
  


 
  
  #######################################################################################
  #######################################################################################
  ############################# 第2种画法#############################
  ####################################################################################### 
  #######################################################################################
  
  ego_list <- list(ego, ego_MF, ego_BP, ego_CC, kegg)
  title_list = c("叠加_GO_BP enrichment barplot", "叠加_GO_CC enrichment barplot", "叠加_ego enrichment barplot", "叠加_GO_MF enrichment barplot", "叠加_KEGG enrichment barplot")
  for (i in 1:length(ego_list)) {
    dt <- as.data.frame(ego_list[[i]])
    # 按照 -log10(pvalue) 排序
    dt <- dt[order(dt$pvalue, decreasing = FALSE), ]
    # 取前 15 条通路
    dt <- head(dt, n = 15)
    #按照通路指定绘图顺序(转换为因子)：
    dt$Description <- factor(dt$Description,
                             levels = rev(dt$Description))
  #在自定义主题中去掉y轴通路标签:
  mytheme2 <- mytheme + theme(axis.text.y = element_blank())
  p1 <- ggplot(data = dt, aes(x = Count, y = Description, fill = -log10(pvalue))) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
    labs(x = "Number of Gene", y = "", title = title_list[i]) +
    geom_text(aes(x = 0.03, #用数值向量控制文本标签起始位置
                  label = Description),
              hjust = 0)+ #hjust = 0,左对齐
    theme_classic() + mytheme2
  p1
  
  
  # 保存为PDF文件
  ggsave(file = paste0(picDir, "/",title_list[i], ".pdf"), plot = p1, width = 6, height = 4.5, units = "in")
  }
  #######################################################################################
  #######################################################################################
  ############################# 第3种画法#############################
  ####################################################################################### 
  #######################################################################################
  
  
  ego_list <- list(ego, ego_MF, ego_BP, ego_CC, kegg)
  title_list = c("第3种_GO_BP enrichment barplot", "第3种_GO_CC enrichment barplot", "第3种_ego enrichment barplot", "第3种_GO_MF enrichment barplot", "第3种_KEGG enrichment barplot")
  for (i in 1:length(ego_list)) {
    dt <- as.data.frame(ego_list[[i]])
    # 按照 -log10(pvalue) 排序
    dt <- dt[order(dt$pvalue, decreasing = FALSE), ]
    # 取前 15 条通路
    dt <- head(dt, n = 15)
    #按照通路指定绘图顺序(转换为因子)：
    dt$Description <- factor(dt$Description,
                             levels = rev(dt$Description))
    #先分别绘制:
    #绘制富集条形图展示-log10p：
    p2 <- ggplot(data = dt,aes(x = Description, y = -log10(pvalue))) +
      geom_bar(stat = "identity", width = 0.8, fill = 'pink') +
      scale_y_continuous(limits = c(0, 25)) +
      labs(x = "", y = "-log10(pvalue)", title = "KEGG enrichment barplot") +
      theme_classic() + mytheme +
      coord_flip()
    p2
    #绘制折线散点图展示count数目：
    p3 <- ggplot(data = dt, aes(x = Description, y = Count, group = 1)) +
      geom_line(color = "grey", cex = 0.5) + #绘制折线
      geom_point(size = 2) + #绘制Count散点
      scale_y_continuous(limits = c(20, 60)) +
      labs(x = "", y = "Count", title = "KEGG enrichment barplot") +
      theme_classic() + mytheme +
      coord_flip()
    p3
    #将不同坐标系图表进行组合：
    p4 <- ggplot() +
      geom_bar(data = dt, aes(x = Description, y = -log10(pvalue)),
               stat = "identity", width = 0.8, fill = 'pink') +
      scale_y_continuous(limits = c(0,25), #主y轴刻度限制
                         breaks = seq(0,25,5), #主y轴刻度线断点
                         expand = c(0,0),
                         sec.axis = sec_axis(~. *2.5, #创建与主轴相对的辅助轴，通过函数计算实现1对1的转换。如主轴范围为0-25，这里我们通过*2.4实现辅助轴范围为0-60
                                             name = 'Count', #辅助轴标题
                                             breaks = seq(0,60,10))) + #辅助轴刻度线断点
      geom_line(data = dt, aes(x = Description, y = Count/2.4, group = 1),
                color = "grey", cex = 0.5) +
      geom_point(data = dt, aes(x = Description, y = Count/2.4),
                 size = 2) +
      labs(x = "", y = "-log10(pvalue)", title = "KEGG enrichment barplot") +
      theme_classic() +
      mytheme +
      coord_flip() #由于geom_line函数绘制折线是按照变量在x轴水平方向顺序进行连接，因此这里把x/y轴对应变量交换了一下，最后再通过坐标翻转回来
    p4
    
  
    # 保存为PDF文件
    ggsave(file = paste0(picDir, "/",title_list[i], ".pdf"), plot = p4, width = 8, height = 6, units = "in")
  }
  
  
  
  
  
  
  # 返回富集分析结果
  return(list(GO_ALL = ego, GO_MF = ego_MF, GO_BP = ego_BP, GO_CC = ego_CC, KEGG = kegg))
}

